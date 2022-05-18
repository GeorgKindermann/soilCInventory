kinterFix:	kinterFix.cc
	g++ kinterFix.cc -O3 -okinterFix -Wall -ffast-math -march=native -funroll-loops -flto -std=c++11 -DRCPP
	strip kinterFix
