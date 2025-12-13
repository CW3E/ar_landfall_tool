## AR Landfall Tool

---

This repository runs calculations and plots for the AR Landfall Tool. Plot types include AR Landfall Contour, AR Landfall Vector, AR Landfall IVT magnitude (control and ensemble-mean) for GEFS, ECMWF, West-WRF, and ECMWF-GEFS. 

Current capabilities includes coastal, foothills, and inland transects for the US-West coast from 25°N to 60°N, coastal, foothills, and inland transects for southern Alaska, and coastal and inland transects for western Alaska.

### To run:

---

To run all three regions with a singularity container:

```bash
## runs plots for GEFS
apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "GEFS" "YYYYMMDDHH"

## runs plots for ECWMF
apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF" "YYYYMMDDHH"

## runs plots for W-WRF
apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "W-WRF" "YYYYMMDDHH"

## runs plots for ECMWF-GEFS
apptainer exec -e --bind /data:/data /data/projects/operations/LandfallTools/ar_landfall_tool/envs/ar_landfall_tool.2025.12.12.sif python /data/projects/operations/LandfallTools/ar_landfall_tool/run_tool.py "ECMWF-GEFS" "YYYYMMDDHH"
```