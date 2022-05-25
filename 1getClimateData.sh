curl --head "https://s3.hub.zamg.ac.at/datahub/index.html?id=5f52dc16-6e25-4a88-9cc1-f7d1e97ab8ff/filelisting&anonymous=true#" -o /dev/null -c /tmp/cookie_zamg.txt

curl -O https://s3.hub.zamg.ac.at/datahub/resources/5f52dc16-6e25-4a88-9cc1-f7d1e97ab8ff/filelisting/SPARTACUS-MONTHLY_DEM_MASK.nc -b /tmp/cookie_zamg.txt --output-dir /dat/meteo/zamg/spartacus/

for YEAR in {1961..2021}
do
    for WHAT in RR Tm
    do
	curl -O https://s3.hub.zamg.ac.at/datahub/resources/5f52dc16-6e25-4a88-9cc1-f7d1e97ab8ff/filelisting/SPARTACUS-MONTHLY_${WHAT}_${YEAR}.nc -b /tmp/cookie_zamg.txt --output-dir /dat/meteo/zamg/spartacus/month
    done
done

