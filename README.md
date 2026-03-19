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
```

---

## Inputs

### Forecast & Hindcast

- NetCDF format  
- Dimensions:
  ```
  time, lat, lon, member
  ```

- Hindcasts may include:
  - **staggered members** (different initialization times per member)
  - **continuous members** (single initialization)

---

### Observations

- NetCDF format  
- Dimensions:
  ```
  time, lat, lon
  ```

---

## Output Structure

### Hindcast & Forecast

```
lead_time, init_date, member, lat, lon
```

### Observations

```
lead_time, init_date, lat, lon
```

---

## Core Concepts

### Lead-Time Representation

All datasets are transformed into a structure where:

- `lead_time = 0` → first day after initialization  
- `init_date` → initialization date of the forecast/hindcast  

This allows consistent analysis across lead times, such as:

- weekly means (week 1, week 2, …)  
- skill vs lead time  
- aggregated forecast periods  

---

### Staggered Hindcasts

The pipeline automatically handles hindcasts where:

- different ensemble members start on different dates  
- missing values (NaNs) mark inactive periods  

Cycle starts are detected using **member-specific NaN transitions**, ensuring correct reconstruction of forecast cycles.

---

### Temporal Aggregation

User-defined aggregation is applied after reshaping:

- `agg_window = 1` → daily  
- `agg_window = 5` → pentads  
- `agg_window = 7` → weekly  

Supported methods:

- `"mean"`
- `"sum"`

---

## Configuration

Optional keyword arguments can be passed to internal steps:

```python
fn.preprocess_forecast(
    ...,
    time_alignment_kwargs=dict(raise_if_missing=False),
    grid_alignment_kwargs=dict(method="bilinear"),
)
```

---

## Logging

Enable logging with:

```python
verbose=True
```

This provides:

- step-by-step processing information  
- clear function boundaries  
- clean error messages  

---

## Notes

- ECMWF precipitation (`tp`) is provided as **running accumulation** and is automatically converted to daily totals  
- Spatial alignment ensures model data are restricted to locations where observations exist  
- Temporal alignment maps observations onto hindcast time structure  
- Designed primarily for **daily data**  

---

## Dependencies

- xarray  
- numpy  
- pandas  
- xesmf  
- netCDF4  

---

## Project Structure

```
functions/
    __init__.py
    functions.py
README.md
```

---

## Future Improvements

- Parallelization for large datasets  
- Support for sub-daily data  
- Built-in verification metrics  
- CLI interface  

---

## Summary

This pipeline solves a non-trivial problem in S2S workflows:

> transforming heterogeneous forecast, hindcast, and observational data  
> into a consistent, lead-time-based structure suitable for analysis

with particular care for:

- staggered hindcasts  
- missing data handling  
- reproducible preprocessing  
