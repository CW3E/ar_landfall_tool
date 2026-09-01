#!/bin/bash

cd /data/projects/website/mirror/htdocs/images/

for LFT in landfalltool*
do
 cd "/data/projects/website/mirror/htdocs/images/${LFT}/v1/ECMWF_ENS/"
 for domain in *
 do
  echo ${LFT} ${domain}
  mv "/data/projects/website/mirror/htdocs/images/${LFT}/v1/ECMWF_ENS/${domain}/2026052800/1/"*"2026052700"* "/data/projects/website/mirror/htdocs/images/${LFT}/v1/ECMWF_ENS/${domain}/2026052700/1/"

 done
done

exit
