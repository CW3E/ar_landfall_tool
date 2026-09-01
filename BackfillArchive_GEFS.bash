#!/bin/bash

start=2026042100
start_epoch=$(date -u -d "${start:0:8} ${start:8:2}:00" +%s)

for i in {0..40}; do
 Date=$(date -u -d "@$((start_epoch - i * 6 * 3600))" +%Y%m%d%H)
 echo "$Date"

# for domain in coast foothills inland intwest
# do
#  for thresh in 150 250 500 750
#  do
#   echo "${domain} ${thresh}"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/${domain}/${Date}/1"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/${domain}/${Date}/1"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__GEFS_50__${domain}__${Date}__1__F384.png"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__GEFS_50__${domain}__${Date}__1__F168.png"
#  done #threshold
#
#  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${Date}/1"
#  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${Date}/1"
#  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivtcontrol__v1__GEFS_50__${domain}__${Date}__1__F384.png"
#  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivtmean__v1__GEFS_50__${domain}__${Date}__1__F384.png"
# done #domain
#
# for domain in intwest
# do
#  for thresh in 100
#  do
#   echo "${domain} ${thresh}"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/${domain}/${Date}/1"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/${domain}/${Date}/1"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__GEFS_50__${domain}__${Date}__1__F384.png"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/US-west/GEFS_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__GEFS_50__${domain}__${Date}__1__F168.png"
#  done #threshold
# done #domain


### AK
 for domain in coast inland
 do
  for thresh in 150 250 500 750
  do
   echo "AK${domain} ${thresh}"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/AK${domain}/${Date}/1"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/AK${domain}/${Date}/1"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/AK/GEFS_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/AK${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__GEFS_50__AK${domain}__${Date}__1__F384.png"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/AK/GEFS_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/AK${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__GEFS_50__AK${domain}__${Date}__1__F168.png"
  done #threshold

  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/AK${domain}/${Date}/1"
  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/AK${domain}/${Date}/1"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/AK/GEFS_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/AK${domain}/${Date}/1/landfalltool_ivtcontrol__v1__GEFS_50__AK${domain}__${Date}__1__F384.png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/AK/GEFS_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/AK${domain}/${Date}/1/landfalltool_ivtmean__v1__GEFS_50__AK${domain}__${Date}__1__F384.png"
 done #domain




### SAK
 for domain in coast inland foothills
 do
  for thresh in 150 250 500 750
  do
   echo "SAK${domain} ${thresh}"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/SAK${domain}/${Date}/1"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/SAK${domain}/${Date}/1"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/SAK/GEFS_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/GEFS_50/SAK${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__GEFS_50__SAK${domain}__${Date}__1__F384.png"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/SAK/GEFS_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/GEFS_50/SAK${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__GEFS_50__SAK${domain}__${Date}__1__F168.png"
  done #threshold

  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/SAK${domain}/${Date}/1"
  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/SAK${domain}/${Date}/1"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/SAK/GEFS_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/GEFS_50/SAK${domain}/${Date}/1/landfalltool_ivtcontrol__v1__GEFS_50__SAK${domain}__${Date}__1__F384.png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/gefs/LandfallTool/${Date}/SAK/GEFS_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/GEFS_50/SAK${domain}/${Date}/1/landfalltool_ivtmean__v1__GEFS_50__SAK${domain}__${Date}__1__F384.png"
 done #domain


done #date
exit
