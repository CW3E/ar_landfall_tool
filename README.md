## AR Landfall Tool

---

This repository runs calculations and plots for the AR Landfall Tool. Plot types include AR Landfall Contour, AR Landfall Vector, AR Landfall IVT magnitude (control and ensemble-mean) for GEFS, ECMWF, West-WRF, and ECMWF-GEFS. 

Current capabilities includes coastal, foothills, and inland transects for the US-West coast from 25°N to 60°N, coastal, foothills, and inland transects for southern Alaska, and coastal and inland transects for western Alaska.

### To run:

---

To run all three regions with a singularity container:

```bash
## runs plots for US-West coast
singularity exec --bind /data:/data,/home/dnash/repos/cw3e_ar-tools:/cw3e_ar-tools,/work:/work,/common:/common -e ar_landfall_tool.sif /opt/conda/bin/python /home/dnash/repos/cw3e_ar-tools/ar_landfall_tool_US-West.py

## runs plots for southern Alaska coast
singularity exec --bind /data:/data,/home/dnash/repos/cw3e_ar-tools:/cw3e_ar-tools,/work:/work,/common:/common -e ar_landfall_tool.sif /opt/conda/bin/python /home/dnash/repos/cw3e_ar-tools/ar_landfall_tool_AK-south.py

## runs plots for western Alaska coast
singularity exec --bind /data:/data,/home/dnash/repos/cw3e_ar-tools:/cw3e_ar-tools,/work:/work,/common:/common -e ar_landfall_tool.sif /opt/conda/bin/python /home/dnash/repos/cw3e_ar-tools/ar_landfall_tool_AK-west.py
```