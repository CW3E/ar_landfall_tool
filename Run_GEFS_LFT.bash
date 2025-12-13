#!/bin/bash

date=`date`
echo "Starting AT "$date

lag=5
yyyy=`date -d '-'$lag' hours' -u +%Y`
mm=`date -d '-'$lag' hours' -u +%m`
dd=`date -d '-'$lag' hours' -u +%d`
hh=`date -d '-'$lag' hours' -u +%H`

mkdir "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/US-west"
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/SAK"
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/AK"


rm -f /data/projects/operations/LandfallTools/figs/US-west/GEFS_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/AK/GEFS_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/SAK/GEFS_LandfallTool*current.png

cd /data/projects/operations/LandfallTools/ar_landfall_tool/


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
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/gefs/v12/LFT/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/AK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/gefs/v12/LFT/AK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/SAK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/gefs/v12/LFT/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/US-west
cp "GEFS_LandfallTool_250_coast_"$yyyy$mm$dd$hh".png" "GEFS_LandfallTool_500_coast_"$yyyy$mm$dd$hh".png" /data/projects/for_sharing/Forecast_Product_Archive/PSing/GEFS/
mv "GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/US-west"
cd /data/projects/operations/LandfallTools/figs/AK
mv "GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/AK"
cd /data/projects/operations/LandfallTools/figs/SAK
mv "GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$yyyy$mm$dd$hh"/SAK"

cd /data/projects/operations/LandfallTools/figs/dProgdT/
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/*.png
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/*.png
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/*.png

X=0
Y=38
while [ $X -lt 169 ];
do
 date=`date -u -d $yyyy$mm$dd' '$hh' -'$X' hours' +%Y%m%d%H`
 for trans in coast foothills inland intwest
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/US-west/GEFS_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/GEFS_LandfallTool_Mean_"$trans"_"$Y".png"
done

 for trans in coast foothills inland
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/SAK/GEFS_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/GEFS_LandfallTool_Mean_"$trans"_"$Y".png"
 done

 for trans in coast inland
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/"$date"/AK/GEFS_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/GEFS_LandfallTool_Mean_"$trans"_"$Y".png"
 done

 let X=X+6
 let Y=Y-1
done

check=1
while [ $check != 0 ]; do
 timeout 220 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/US-west/ /data/projects/website/mirror/htdocs/images/gefs/v12/dprog/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

check=1
while [ $check != 0 ]; do
 timeout 220 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/AK/ /data/projects/website/mirror/htdocs/images/gefs/v12/dprog/AK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

check=1
while [ $check != 0 ]; do
 timeout 220 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/GEFS/SAK/ /data/projects/website/mirror/htdocs/images/gefs/v12/dprog/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

date=`date`
echo "Finished at "$date


exit





