date=$1
hh=$2
URLp="https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?file=gfs.t"${hh}"z.pgrb2.0p25.f168&lev_surface=on&var_APCP=on&leftlon=-175&rightlon=-125&toplat=75&bottomlat=40&dir=%2Fgfs."${date}"%2F"${hh}"%2Fatmos"
curl "$URLp" -o precip.grb