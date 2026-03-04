#!/bin/bash

date=$1
hh=$2

echo "DOWNLOADING GFS PRECIP DATA FOR ${date} ${hh}"


#### OpenDAP retirement 2026 02 ####
#URLp="https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t"${hh}"z.pgrb2.0p25.f168&lev_surface=on&var_APCP=on&leftlon=-175&rightlon=-125&toplat=75&bottomlat=40&dir=%2Fgfs."${date}"%2F"${hh}"%2Fatmos"
URLp="https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.${date}%2F${hh}%2Fatmos&file=gfs.t${hh}z.pgrb2.0p25.f168&var_APCP=on&lev_surface=on&subregion=&toplat=75&leftlon=-175&rightlon=-125&bottomlat=40"

echo "DOWNLOADING ${URLp}"

curl "$URLp" -o precip_GFS.grb
