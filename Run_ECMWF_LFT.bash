#!/bin/bash

date=`date`
echo "Starting AT "$date

lag=7
yyyy=`date -d '-'$lag' hours' -u +%Y`
mm=`date -d '-'$lag' hours' -u +%m`
dd=`date -d '-'$lag' hours' -u +%d`
hh=`date -d '-'$lag' hours' -u +%H`

mkdir "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/US-west/"
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/SAK/"
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/AK"

rm -f /data/projects/operations/LandfallTools/figs/US-west/ECMWF_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/AK/ECMWF_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/SAK/ECMWF_LandfallTool*current.png

cd /data/projects/operations/LandfallTools/ar_landfall_tool/


filename="/data/projects/derived_products/ECMWF_IVT/Ensemble/IVT_EC_"$yyyy$mm$dd$hh".nc"
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

singularity exec --bind /data:/data -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF" "$yyyy$mm$dd$hh"

cd /data/projects/operations/LandfallTools/figs/US-west
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/AK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/AK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/SAK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

echo "/data/projects/operations/LandfallTools/figs/US-west/"
ls /data/projects/operations/LandfallTools/figs/US-west/
echo "/data/projects/operations/LandfallTools/figs/AK/"
ls /data/projects/operations/LandfallTools/figs/AK/
echo "/data/projects/operations/LandfallTools/figs/SAK/"
ls /data/projects/operations/LandfallTools/figs/SAK/


cp "/data/projects/operations/LandfallTools/figs/US-west/ECMWF_LandfallTool_250_coast_"$yyyy$mm$dd$hh".png" "ECMWF_LandfallTool_500_coast_"$yyyy$mm$dd$hh".png" /data/projects/for_sharing/Forecast_Product_Archive/PSing/ECMWF/
mv "/data/projects/operations/LandfallTools/figs/US-west/ECMWF_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/US-west/"
mv "/data/projects/operations/LandfallTools/figs/AK/ECMWF_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/AK/"
mv "/data/projects/operations/LandfallTools/figs/SAK/ECMWF_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/SAK/"

cd /data/projects/operations/LandfallTools/figs/dProgdT/
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/*.png
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/*.png
rm -f /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/*.png

X=0
Y=24
while [ $X -lt 169 ];
do
 date=`date -u -d $yyyy$mm$dd' '$hh' -'$X' hours' +%Y%m%d%H`
 for trans in coast foothills inland intwest
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/US-west/ECMWF_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ECMWF_LandfallTool_Mean_"$trans"_"$Y".png"
done

 for trans in coast foothills inland
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/SAK/ECMWF_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ECMWF_LandfallTool_Mean_"$trans"_"$Y".png"
 done

 for trans in coast inland
 do
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_150_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_150_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_250_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_250_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_500_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_500_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_750_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_750_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_control_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_Control_"$trans"_"$Y".png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$date"/AK/ECMWF_LandfallTool_ensemble_mean_"$trans"_"$date".png" "/data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ECMWF_LandfallTool_Mean_"$trans"_"$Y".png"

 done

 let X=X+12
 let Y=Y-1
done

check=1
while [ $check != 0 ]; do
 timeout 120 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/US-west/ /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/dprog/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

check=1
while [ $check != 0 ]; do
 timeout 120 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/AK/ /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/dprog/AK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

check=1
while [ $check != 0 ]; do
 timeout 120 rsync --delete-before -avih /data/projects/operations/LandfallTools/figs/dProgdT/ECMWF/SAK/ /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/dprog/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

date=`date`
echo "Finished at "$date


exit





