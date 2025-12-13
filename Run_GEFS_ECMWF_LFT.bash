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

rm -f /data/projects/operations/LandfallTools/figs/US-west/ECMWF-GEFS_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/AK/ECMWF-GEFS_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/SAK/ECMWF-GEFS_LandfallTool*current.png

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

apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF-GEFS" "$yyyy$mm$dd$hh"


cd /data/projects/operations/LandfallTools/figs/US-west
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF-GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/AK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF-GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/AK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/SAK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih ECMWF-GEFS_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/ECMWF/ensemble/LandfallTool/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done


cd /data/projects/operations/LandfallTools/figs/US-west
mv "ECMWF-GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/US-west/"
cd /data/projects/operations/LandfallTools/figs/AK
mv "ECMWF-GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/AK/"
cd /data/projects/operations/LandfallTools/figs/SAK
mv "ECMWF-GEFS_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/"$yyyy$mm$dd$hh"/SAK/"


date=`date`
echo "Finished at "$date


exit





