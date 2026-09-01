#!/bin/bash

start=2026042012
start_epoch=$(date -u -d "${start:0:8} ${start:8:2}:00" +%s)

for i in {0..39}; do
 Date=$(date -u -d "@$((start_epoch - i * 12 * 3600))" +%Y%m%d%H)

 echo "$Date"

# for domain in coast foothills inland intwest
# do
#  for thresh in 150 250 500 750
#  do
#   echo "${domain} ${thresh}"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/${domain}/${Date}/1"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/${domain}/${Date}/1"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__ECMWF_ENS__${domain}__${Date}__1__F384.png"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__ECMWF_ENS__${domain}__${Date}__1__F168.png"
#  done #threshold
#
#  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/${domain}/${Date}/1"
#  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/${domain}/${Date}/1"
#  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivtcontrol__v1__ECMWF_ENS__${domain}__${Date}__1__F384.png"
#  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivtmean__v1__ECMWF_ENS__${domain}__${Date}__1__F384.png"
# done #domain

# for domain in intwest
# do
#  for thresh in 100
#  do
#   echo "${domain} ${thresh}"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/${domain}/${Date}/1"
#   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/${domain}/${Date}/1"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__ECMWF_ENS__${domain}__${Date}__1__F384.png"
#   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/US-west/ECMWF_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__ECMWF_ENS__${domain}__${Date}__1__F168.png"
#  done #threshold
# done #domain


### AK
 for domain in coast inland
 do
  for thresh in 150 250 500 750
  do
   echo "AK${domain} ${thresh}"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/AK${domain}/${Date}/1"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/AK${domain}/${Date}/1"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/AK/ECMWF_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/AK${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__ECMWF_ENS__AK${domain}__${Date}__1__F384.png"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/AK/ECMWF_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/AK${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__ECMWF_ENS__AK${domain}__${Date}__1__F168.png"
  done #threshold

  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/AK${domain}/${Date}/1"
  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/AK${domain}/${Date}/1"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/AK/ECMWF_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/AK${domain}/${Date}/1/landfalltool_ivtcontrol__v1__ECMWF_ENS__AK${domain}__${Date}__1__F384.png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/AK/ECMWF_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/AK${domain}/${Date}/1/landfalltool_ivtmean__v1__ECMWF_ENS__AK${domain}__${Date}__1__F384.png"
 done #domain

### SAK
 for domain in coast inland foothills
 do
  for thresh in 150 250 500 750
  do
   echo "SAK${domain} ${thresh}"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/SAK${domain}/${Date}/1"
   mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/SAK${domain}/${Date}/1"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/SAK/ECMWF_LandfallTool_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_probability/v1/ECMWF_ENS/SAK${domain}/${Date}/1/landfalltool_ivt${thresh}_probability__v1__ECMWF_ENS__SAK${domain}__${Date}__1__F384.png"
   cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/SAK/ECMWF_LandfallTool_Vectors_${thresh}_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivt${thresh}_vectors/v1/ECMWF_ENS/SAK${domain}/${Date}/1/landfalltool_ivt${thresh}_vectors__v1__ECMWF_ENS__SAK${domain}__${Date}__1__F168.png"
  done #threshold

  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/SAK${domain}/${Date}/1"
  mkdir -p -m 2775 "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/SAK${domain}/${Date}/1"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/SAK/ECMWF_LandfallTool_control_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtcontrol/v1/ECMWF_ENS/SAK${domain}/${Date}/1/landfalltool_ivtcontrol__v1__ECMWF_ENS__SAK${domain}__${Date}__1__F384.png"
  cp "/data/downloaded/Forecasts/ARPortal_Archive/ECMWF/ensemble/landfall_images/${Date}/SAK/ECMWF_LandfallTool_ensemble_mean_${domain}_${Date}.png" "/data/projects/website/mirror/htdocs/images/landfalltool_ivtmean/v1/ECMWF_ENS/SAK${domain}/${Date}/1/landfalltool_ivtmean__v1__ECMWF_ENS__SAK${domain}__${Date}__1__F384.png"
 done #domain


done #date
exit

