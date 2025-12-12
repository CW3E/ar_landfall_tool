## AR Landfall Tool

---

This repository runs calculations and plots for the AR Landfall Tool. Plot types include AR Landfall Contour, AR Landfall Vector, AR Landfall IVT magnitude (control and ensemble-mean) for GEFS, ECMWF, West-WRF, and ECMWF-GEFS. 

Current capabilities includes coastal, foothills, and inland transects for the US-West coast from 25°N to 60°N, coastal, foothills, and inland transects for southern Alaska, and coastal and inland transects for western Alaska.

### To run:

---

To run all three regions with a singularity container:

```bash
## runs plots for GEFS
singularity exec --bind /data:/data,/home:/home,/work:/work,/common:/common -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "GEFS" "YYYYMMDDHH"

## runs plots for ECWMF
singularity exec --bind /data:/data,/home:/home,/work:/work,/common:/common -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF" "YYYYMMDDHH"

## runs plots for W-WRF
singularity exec --bind /data:/data,/home:/home,/work:/work,/common:/common -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "W-WRF" "YYYYMMDDHH"

## runs plots for ECMWF-GEFS
singularity exec --bind /data:/data,/home:/home,/work:/work,/common:/common -e /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.sif /opt/conda/bin/python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF-GEFS" "YYYYMMDDHH"
```