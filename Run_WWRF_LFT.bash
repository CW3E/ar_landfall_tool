#!/bin/bash

date=`date`
echo "Starting AT "$date

lag=13
yyyy=`date -d '-'$lag' hours' -u +%Y`
mm=`date -d '-'$lag' hours' -u +%m`
dd=`date -d '-'$lag' hours' -u +%d`
hh=`date -d '-'$lag' hours' -u +%H`

mkdir "/data/downloaded/Forecasts/ARPortal_Archive/WWRF/ensemble/"$yyyy$mm$dd$hh
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/WWRF/ensemble/"$yyyy$mm$dd$hh"/US-west"
mkdir "/data/downloaded/Forecasts/ARPortal_Archive/WWRF/ensemble/"$yyyy$mm$dd$hh"/SAK"



rm -f /data/projects/operations/LandfallTools/figs/US-west/W-WRF_LandfallTool*current.png
rm -f /data/projects/operations/LandfallTools/figs/SAK/W-WRF_LandfallTool*current.png


cd /data/projects/operations/LandfallTools/ar_landfall_tool/


filename="/data/downloaded/WWRF-NRT/2025-2026/Ensemble_IVT/IVT_WWRF_"$yyyy$mm$dd$hh".nc"
last_mtime=""

while true; do
  if [[ -e "$filename" ]]; then
    current_mtime=$(stat --format="%Y" "$filename")

    if [[ "$last_mtime" != "$current_mtime" ]]; then
      echo "$(date): Detected new or updated $filename. Checking size stability..."

      while true; do
        filesize1=$(stat --format="%s" "$filename")
        sleep 5
        filesize2=$(stat --format="%s" "$filename")

        if [[ "$filesize1" == "$filesize2" ]]; then
          echo "$(date): $filename ready for processing"

###### Processing Code ######
date=`date`
echo "STARTING MAKING LFTs at "$date

singularity exec --bind /data:/data -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "W-WRF" "$yyyy$mm$dd$hh"

cd /data/projects/operations/LandfallTools/figs/US-west
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih W-WRF_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/wwrf/images/ensemble/LFT/US-west/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/SAK/
check=1
while [ $check != 0 ]; do
timeout 120 rsync --ignore-missing-args -avih W-WRF_LandfallTool*current.png /data/projects/website/mirror/htdocs/images/wwrf/images/ensemble/LFT/SAK/
 check=$?
 if [ $check != 0 ]; then
  echo "Transfer failed"
  sleep 10
 fi
done

cd /data/projects/operations/LandfallTools/figs/US-west
mv "W-WRF_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/WWRF/ensemble/"$yyyy$mm$dd$hh"/US-west"
cd /data/projects/operations/LandfallTools/figs/AK
mv "W-WRF_LandfallTool"*$yyyy$mm$dd$hh".png" "/data/downloaded/Forecasts/ARPortal_Archive/WWRF/ensemble/"$yyyy$mm$dd$hh"/SAK"

date=`date`
echo "Finished at "$date
#############################
        last_mtime="$current_mtime"
        break
        else
          echo "$(date): $filename is still being written (file size is changing)."
          sleep 30
        fi
      done
    else
      echo "$(date): $filename has not been updated since last processing."
    fi
  else
    echo "$(date): $filename does not exist."
  fi

 sleep 60
done


exit





