#!/bin/bash

date=`date`
echo "Starting AT "$date

lag=5
yyyy=`date -d '-'$lag' hours' -u +%Y`
mm=`date -d '-'$lag' hours' -u +%m`
dd=`date -d '-'$lag' hours' -u +%d`
hh=`date -d '-'$lag' hours' -u +%H`

cd /data/projects/operations/LandfallTools/ar_landfall_tool/
unalias cp

filename="/data/projects/derived_products/GEFS_IVT/data/GEFS_IVT_"$yyyy$mm$dd$hh".nc"
while true; do
  if [[ -e "$filename" ]]; then                   # Check if the file exists
    filesize1=$(stat --format="%s" "$filename")   # Get the file size
    sleep 5                                      # Wait a few seconds
    filesize2=$(stat --format="%s" "$filename")   # Get the file size again

    if [[ "$filesize1" == "$filesize2" ]]; then  # Compare file sizes
      break                                      # Exit the loop if file size is not changing
    else
      echo $filename" is still being written (file size is changing)."
      sleep 30
    fi
  else
    echo $filename" does not exist."
    sleep 30
  fi
done

echo $filename" ready for processing"

date=`date`
echo "STARTING MAKING LFTs at "$date

apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "GEFS" "$yyyy$mm$dd$hh"

cd /data/projects/operations/LandfallTools/figs/US-west
chmod 664 landfalltoolpng
for domain in coast foothills inland intwest
do
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mv "landfalltool_ivt150_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt150_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtmean__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtcontrol__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
done
wait

for domain in intwest 
do
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mv "landfalltool_ivt100_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt100_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
done
wait

cd /data/projects/operations/LandfallTools/figs/SAK
chmod 664 landfalltoolpng
for domain in SAKcoast SAKfoothills SAKinland
do
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mv "landfalltool_ivt150_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt150_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtmean__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtcontrol__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
done
wait

cd /data/projects/operations/LandfallTools/figs/AK
chmod 664 landfalltoolpng
for domain in AKcoast AKinland
do
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1"
 mv "landfalltool_ivt150_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt150_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt250_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt500_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivt750_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtmean__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
 mv "landfalltool_ivtcontrol__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}"* "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1" &
done
wait

cp /data/projects/operations/LandfallTools/figs/US-west/GEFS_LandfallTool_250_coast_current.png "/data/projects/for_sharing/Forecast_Product_Archive/PSing/GEFS/GEFS_LandfallTool_250_coast_${yyyy}${mm}${dd}${hh}.png" 
cp /data/projects/operations/LandfallTools/figs/US-west/GEFS_LandfallTool_500_coast_current.png "/data/projects/for_sharing/Forecast_Product_Archive/PSing/GEFS/GEFS_LandfallTool_500_coast_${yyyy}${mm}${dd}${hh}.png"

for domain in coast foothills inland intwest SAKcoast SAKfoothills SAKinland AKcoast AKinland
do
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt150_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_probability/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt150_probability__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt250_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_probability/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt250_probability__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt500_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_probability/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt500_probability__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt750_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_probability/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt750_probability__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt150_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt150_vectors/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt150_vectors__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt250_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt250_vectors/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt250_vectors__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt500_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt500_vectors/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt500_vectors__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt750_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt750_vectors/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt750_vectors__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivtmean__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/latest/1/landfalltool_ivtmean__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivtcontrol__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/latest/1/landfalltool_ivtcontrol__v1__GEFS_50__${domain}__latest__1_ F384.png" &
 wait
done

for domain in intwest
do
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_probability/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt100_probability__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F384.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_probability/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt100_probability__v1__GEFS_50__${domain}__latest__1__F384.png" &
 cp "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_vectors/v1/GEFS_50/${domain}/${yyyy}${mm}${dd}${hh}/1/landfalltool_ivt100_vectors__v1__GEFS_50__${domain}__${yyyy}${mm}${dd}${hh}__1__F168.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt100_vectors/v1/GEFS_50/${domain}/latest/1/landfalltool_ivt100_vectors__v1__GEFS_50__${domain}__latest__1__F168.png" &
done
wait



date=`date`
echo "Finished at "$date


exit





