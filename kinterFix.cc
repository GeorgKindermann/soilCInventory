#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <unordered_map>

using namespace std;

#ifndef RCPP
  #include <Rcpp.h>
  #include <cstring>
  using namespace Rcpp;
#endif

static void usage (char argv[]) {
  cout << argv << " - Klima Interpolation, (C) 2022 Georg Kindermann" << endl;
  cout << "Usage: " << argv << " <OPTIONS>" << endl;
  cout << "Options:" << endl;
  cout << "  -o, --output NAME         send output to NAME instead of out.txt" << endl;
  cout << "  -m, --measure NAME        read measurements from NAME instead of mdata.txt" << endl;
  cout << "  -p, --position NAME       read point positions from NAME instead of pdata.txt" << endl;
  cout << "  -s, --stations n          Minimum number of stations to use (def=8)" << endl;
  cout << "  -r, --radius x            Maximum distance of stations to use (def=1.)" << endl;
  cout << "  -fx x                     Factor to multiply x-distances (def=1.)" << endl;
  cout << "  -fy x                     Factor to multiply y-distances (def=1.)" << endl;
  cout << "  -fz x                     Factor to multiply z-distances (def=1.)" << endl;
  cout << "  -gd x                     Default value of gradient (def=0.)" << endl;
  cout << "  -gl x                     Lowest value of gradient (def=-INF)" << endl;
  cout << "  -gh x                     Highest value of gradient (def=INF)" << endl;
  cout << "  -rl x                     Lowest value of result (def=-INF)" << endl;
  cout << "  -rh x                     Highest value of result (def=-INF)" << endl;
  cout << "  -h, --help                display this help and exit" << endl;
  cout << "  -V, --version             output version information and exit" << endl;
  cout << "Formats:" << endl;
  cout << "  Output: x y z date interpolatedValue <...>" << endl;
  cout << "  Measurement: x y z date measuredValue" << endl;
  cout << "    Measurementfile needs to be date sorted" << endl;
  cout << "  Point Position: x y z <...>" << endl;
  exit(0);
}

struct point {
  double x;  //Position
  double y;
  double z;
  double v;  //Value Messwert oder Interpolation
};

struct posD {
  double x;  //Position
  double y;
  double z;
  bool operator==(const posD &other) const
  {return (x == other.x && y == other.y && z == other.z);}
  bool operator<(const posD &other) const {
    return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);}
};

struct posI {
  int x;  //Position
  int y;
  int z;
  bool operator==(const posI &other) const {
    return (x == other.x && y == other.y && z == other.z);}
  bool operator<(const posI &other) const {
    return std::tie(x, y, z) < std::tie(other.x, other.y, other.z);}
};

namespace std {
  template <>
  struct hash<posD> {
    std::size_t operator()(const posD& k) const
    {
      using std::size_t;
      using std::hash;
      return ((hash<double>()(k.x)
               ^ (hash<double>()(k.y) << 1)) >> 1)
	^ (hash<double>()(k.z) << 1);
    }
  };

  template <>
  struct hash<posI> {
    std::size_t operator()(const posI& k) const
    {
      using std::size_t;
      using std::hash;
      return ((hash<int>()(k.x)
               ^ (hash<int>()(k.y) << 1)) >> 1)
	^ (hash<int>()(k.z) << 1);
    }
  };
}

bool cmp(double *a, double *b ) {
   return(*a < *b);
} 

int interpol(//Ausgesuchte naheliegende Stationen
	     vector<vector<map<posD, double>::iterator> > selectedStations
	     , vector<vector<double> > selectedStationsDist
	     , vector<point> &est//est..Points where Values sholud be estimated
	     , double fx = 1.  //Multiplikator fuer X-Entfernungen
	     , double fy = 1.  //Multiplikator fuer Y-Entfernungen
	     , double fz = 1.  //Multiplikator fuer Z-Entfernungen
	     , double slDef = 0. //default Slope
	     , double slMin = -INFINITY //minimale Slope
	     , double slMax = INFINITY  //maximale Slope
	     , double min = -INFINITY   //minimum der Schaetzung
	     , double max = INFINITY   //maximum der Schaetzung
	     ) {
  if(selectedStations.size() == 0) {
    for(vector<point>::iterator pe=est.begin(); pe!=est.end(); ++pe) {
      pe->v = NAN;
    }
  } else {
    for(unsigned int i=0; i<est.size(); ++i) {
      double r = NAN;
      if(selectedStations[i].size() > 0) {
	//Beobachtungen genau auf dem gesuchten Punkt
	if(selectedStationsDist[i][0] == 0) {
	  r = selectedStations[i][0]->second;
	} else {
	  if(selectedStationsDist[i].size() == 1){
	    r = selectedStations[i][0]->second + slDef*(est[i].z - selectedStations[i][0]->first.z);
	  } else {
	    unsigned int n = selectedStationsDist[i].size();
	    double maxdist = sqrt(selectedStationsDist[i][n-1]*selectedStationsDist[i][n-1] * (n+0.5) / double(n));
	    maxdist = 1./maxdist;
	    double h=0.; double x=0.; double g=0.;
	    for(unsigned int j=0; j<n; ++j) {
	      h += selectedStations[i][j]->first.z * (1./selectedStationsDist[i][j] - maxdist);
	      x += selectedStations[i][j]->second * (1./selectedStationsDist[i][j] - maxdist);
	      g += (1./selectedStationsDist[i][j] - maxdist);
	    }
	    h /= g; x /= g;
	    double sz=0.; double sn=0.;
	    {
	      for(unsigned int j=0; j<n; ++j) {
		sz *= (selectedStations[i][j]->first.z - h) * (selectedStations[i][j]->second - x) * (1./selectedStationsDist[i][j] - maxdist);
		sn += (selectedStations[i][j]->first.z - h) * (selectedStations[i][j]->first.z - h) * (1./selectedStationsDist[i][j] - maxdist);
	      }
	    }
	    double b = sz/sn;
	    if(b != b) {b = slDef;}
	    if(b < slMin) {b = slMin;}
	    if(b > slMax) {b = slMax;}
	    double a = x - b*h;
	    r = a + b*est[i].z;
	  }
	}
	if(r < min) {r = min;}
	if(r > max) {r = max;}
      }
      est[i].v = r;
      //est[i].v = selectedStations[i][0]->second; 
    }
  }
  return(0);
}

int mainX(int argc, char *argv[]) {
  string output = "out.txt";
  string measure = "mdata.txt";
  string position = "pdata.txt";
  int stations = 8;
  double radius = 1.;
  double fx = 1.;
  double fy = 1.;
  double fz = 1.;
  double gd = 0.;
  double gl = -INFINITY;
  double gh = INFINITY;
  double rl = -INFINITY;
  double rh = INFINITY;
  for(int i=1;i<argc;++i) {
    if(argv[i][0] == '-') {
      string arg = argv[i];
      if(arg == "-o" || arg == "--output") {++i; output = argv[i];}
      else if(arg == "-m" || arg == "--measure") {++i; measure = argv[i];}
      else if(arg == "-p" || arg == "--position") {++i; position = argv[i];}
      else if(arg=="-s" || arg == "--stations") {++i; stations = atoi(argv[i]);}
      else if(arg=="-r" || arg == "--radius") {++i; radius = atof(argv[i]);}
      else if(arg=="-fx") {++i; fx = atof(argv[i]);}
      else if(arg=="-fy") {++i; fy = atof(argv[i]);}
      else if(arg=="-fz") {++i; fz = atof(argv[i]);}
      else if(arg=="-gd") {++i; gd = atof(argv[i]);}
      else if(arg=="-gl") {++i; gl = atof(argv[i]);}
      else if(arg=="-gh") {++i; gh = atof(argv[i]);}
      else if(arg=="-rl") {++i; rl = atof(argv[i]);}
      else if(arg=="-rh") {++i; rh = atof(argv[i]);}
      if(arg == "-h" || arg == "--help") {usage(argv[0]);}
      if(arg == "-V" || arg == "--version") {cout << argv[0] << " Version 0.1a, (C) 2022 Georg Kindermann" << endl; exit(0);}
    }
  }
  if(stations < 1) {cerr << "To few Stations" << endl; usage(argv[0]);}
  if(radius <= 0.) {cerr << "To small radius" << endl; usage(argv[0]);}
  if(gh < gl) {double t=gl; gl=gh; gh=t;}
  if(gd < gl) {gd = gl;}
  if(gd > gh) {gd = gh;}
  if(rh < rl) {double t=rl; rl=rh; rh=t;}

  //Punkte zu denen hininterpoliert werden soll einlesen
  vector<point> obs;
  vector<point> est;
  vector<string> beg;
  vector<string> tail;
  {
    ifstream finp(position.c_str());
    if(!finp) {cerr << "Can not open " << position << endl; usage(argv[0]);}
    point tmp;
    string line;
    while(getline(finp, line)) {
      if(line[0] != '#') {
	stringstream ss(line);
	double d[3];
	string store;
	for(int i=0; i<3; ++i) {
	  string tmp;
	  ss >> tmp;
	  if(i>0) {store += "\t";}
	  store += tmp;
	  istringstream inStream(tmp);
	  inStream >> d[i];
	}
	beg.push_back(store);
	{ char drop;
	  ss.get(drop);
	}
	string rest;
	getline(ss, rest);
	tail.push_back(rest);
	tmp.x = d[0];
	tmp.y = d[1];
	tmp.z = d[2];
	tmp.v = NAN;
	est.push_back(tmp);
      }
    }
    finp.close();
  }

  //Alle Positionen der Stationen einlesen
  map<posD, double> stationswerte;
  //unordered_map<posI, double> stationswerte;
  {
    ifstream finm(measure.c_str());
    if(!finm) {cerr << "Can not open " << measure << endl; usage(argv[0]);}
    string line;
    while(finm.good()) {
      getline(finm, line);
      if(line[0] != '#') {
	stringstream ss(line);
	posD pos;
	ss >> pos.x; ss >> pos.y; ss >> pos.z;
	stationswerte.insert({pos, NAN});
      }
    }
  }
  
  //Naeheste Nachbarn suchen und merken
  vector<vector<map<posD, double>::iterator> > selectedStations;
  vector<vector<double> > selectedStationsDist;
  {
    double *dist = new double[stationswerte.size()];
    double **pdist = new double*[stationswerte.size()];
    map<posD, double>::iterator *stIt = new map<posD, double>::iterator[stationswerte.size()];
    for(vector<point>::iterator pe=est.begin(); pe!=est.end(); ++pe) {
      int k=0;
      for(map<posD, double>::iterator po=stationswerte.begin(); po!=stationswerte.end(); ++po){
	dist[k] = sqrt((pe->x - po->first.x)*(pe->x - po->first.x)*fx*fx
		       + (pe->y - po->first.y)*(pe->y - po->first.y)*fy*fy
		       + (pe->z - po->first.z)*(pe->z - po->first.z)*fz*fz);
	pdist[k] = &dist[k];
	stIt[k] = po;
	++k;
      }
      double radMin = radius;
      unsigned int statMin = stations;
      //sort(pdist, pdist + stationswerte.size(), cmp);
      unsigned int middle = statMin;
      if(middle > stationswerte.size()) {middle = stationswerte.size();}
      partial_sort(pdist, pdist + middle, pdist + stationswerte.size(), cmp);
      vector<unsigned int> index;  //Find used Stations
      vector<map<posD, double>::iterator> selStations;
      vector<double> selStationsDist;
      for(unsigned int station=0; station<middle; ++station) {
	unsigned int idx = pdist[station] - dist;
	//Es gibt Beobachtungen genau auf dem gesuchten Punkt -> break
	if(index.size() > 0 && dist[index[0]] == 0. && dist[idx] > 0) {break;}
	if(dist[idx] < radMin && index.size()<statMin) {
	  index.push_back(idx);
	  selStations.push_back(stIt[idx]);
	  selStationsDist.push_back(dist[idx]);
	} else {break;}
      }
      selectedStations.push_back(selStations);
      selectedStationsDist.push_back(selStationsDist);
    }
    delete[] dist;
    delete[] pdist;
    delete[] stIt;
  }
  
  if (est.size() > 0) {
    ifstream finm(measure.c_str());
    if(!finm) {cerr << "Can not open " << measure << endl; usage(argv[0]);}
    ofstream fout(output.c_str(), ios_base::app);
    if(!fout) {cerr << "Can not open " << output << endl; usage(argv[0]);}

    string line;
    bool first = true;
    string cdate;
    string date;
    posD pos;
    double val;
    while(finm.good()) {
      getline(finm, line);
      if(line[0] != '#') {
	stringstream ss(line);
	ss >> pos.x; ss >> pos.y; ss >> pos.z;
	ss >> date;
	ss >> val;
	if(first == true) {first = false; cdate = date;}
	if(cdate != date || !finm.good()) {
	  interpol(selectedStations, selectedStationsDist, est, fx, fy, fz, gd, gl, gh, rl, rh);
	  for(unsigned int i=0; i<est.size(); ++i) {
	    fout << /*beg[i] << "\t" <<*/ cdate << "\t" << est[i].v;
	    if(tail[i].size() > 0) {fout << "\t" << tail[i];}
	    fout << endl;
	  }
	  cdate = date;
	  for(map<posD, double>::iterator po=stationswerte.begin(); po!=stationswerte.end(); ++po) {po->second = NAN;}
	}
	auto search = stationswerte.find(pos);
	if(search != stationswerte.end()) {search->second = val;
	} else {cerr << "Station nicht gefunden" << endl;}
      }
    }
    finm.close();
  }
  return(0); 
}

#ifdef RCPP
int main(int argc, char *argv[]) {
  return mainX(argc, argv);
}
#endif

#ifndef RCPP
// [[Rcpp::export]]
void kinterFix(string arg) {
  vector<string> argvX;
  {
    stringstream ss(arg);
    string tmp;
    while(getline(ss, tmp, ' ')) argvX.push_back(tmp);
  }
  int argc = argvX.size();
  char** argv = new char*[argc]; 
  for(int i = 0; i < argc; ++i) {
    char* s_cstr = new char[argvX[i].size()+1];
    strcpy(s_cstr, argvX[i].c_str());
    argv[i] = s_cstr;
  }
  mainX(argc, argv); 
}
#endif

