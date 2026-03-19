# S2S Forecast Preprocessing Pipeline

A Python pipeline for preprocessing **forecast, hindcast, and observational datasets** into a unified **lead-time framework** for subseasonal-to-seasonal (S2S) analysis and verification.

---

## Overview

This pipeline performs end-to-end preprocessing required to make forecast, hindcast, and observational datasets directly comparable.

The workflow includes:

1. Reading forecast, hindcast, and observation NetCDF files  
2. Converting accumulated variables to daily totals (when required)  
3. Spatial grid alignment between model data and observations  
4. Temporal alignment of observations with hindcast structure  
5. Reshaping datasets into a **lead-time oriented format**  
6. Temporal aggregation over user-defined windows (e.g. 5-day, 7-day)

The final outputs share a **consistent structure**, enabling straightforward comparison across datasets.

---

## Key Features

- Works with both **staggered and non-staggered hindcasts**
- Handles **member-specific initialization times**
- Robust handling of **NaNs and cycle boundaries**
- Converts accumulated variables (e.g. ECMWF precipitation) to daily values
- Produces data ready for:
  - forecast verification  
  - skill analysis  
  - lead-time dependent diagnostics  

---

## Usage

```python
import functions as fn

hc, fc, obs = fn.preprocess_forecast(
    nominal_date="2020-01-01",
    download_dir="/data",
    target_domain="global",
    fcst_model="ECMWF",
    fcst_var="tp",
    agg_window=7,
    agg_method="mean",
    obs_file="/data/obs.nc",
    obs_var="precip",
    verbose=True
)
