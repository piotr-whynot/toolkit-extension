# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# download sub-seasonal forecast data from S2S prediction project
from acacia_s2s_toolkit import argument_check, argument_output, ecdsAPI_requests, download_S2Stc_tracks, download_S2SMJO_indices, grid_functions
import os
import sys
import datetime
import numpy as np

# this is for sphinx - only functions listed here will have entries in readthedocs API
__all__ = ["download_forecast"]

class SuppressOutput:
    """Context manager to silence stdout/stderr (for ECMWF WebAPI logs)."""
    def __enter__(self):
        self._stdout = sys.stdout
        self._stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')
        return self

    def __exit__(self, exc_type, exc, tb):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self._stdout
        sys.stderr = self._stderr
        # If an exception happened, return False so Python re-raises it (with logs visible next time)
        return False

def download_forecast(variable,
                      model='ECMWF',
                      fcdate=None,
                      leadtime_hour=None,
                      start_lt=0,
                      end_lt=None,
                      period='1D',
                      country_name=None,
                      region_name=None,
                      bbox_bounds=[90, -180, -90, 180],
                      plevs=None,
                      filename=None,
                      data_save_dir=None,
                      data_format="netcdf",
                      grid="1.5x1.5",
                      fc_enslags=None,
                      overwrite=False,
                      verbose=True,
                      cleanup=True):

    """
    Overarching function that will download forecast data from ECDS or aux.S2S

    Parameters
    ----------
    model : str
        Forecast model name.
    variable : str
        Variable to download.
    model : str
        Forecast model name. Default = ECMWF
    fcdate : str, optional
        Forecast initialization date in YYYYMMDD format. If None, today's UTC
        date is used and rolled back until a valid forecast date is found.
    plevs : int or list[int], optional
        Pressure level(s), where applicable.
    country_name : str, optional
        Country name used to derive bounds automatically from Natural Earth
        country polygons. If provided, this overrides both region_name and
        bbox_bounds.
    region_name : str, optional
        Predefined domain name. If provided, this overrides bbox_bounds unless
        country_name is also provided.
    bbox_bounds : list[float]
        Bounding box in the form [north, west, south, east].
    filename : str, optional
        Output filename without or with extension.
    data_save_dir : str, optional
        Directory where the downloaded file will be saved.
    data_format : str
        Output format, e.g. "netcdf" or "grib".
    grid : str
        Grid resolution as "dlonxdlat", e.g. "1.5x1.5".
    leadtime_hour : optional
        Forecast lead hours passed to downstream request logic.
    fc_enslags : optional
        Forecast ensemble lag selection passed to downstream request logic.
    overwrite : bool
        If False and file exists already, skip the download.
    verbose : bool
        If True, print diagnostic information.

    Returns
    -------
    str
        Path to the saved file.
    """
    
    # Normalise grid for downstream request code
    grid_for_request = grid.replace("x", "/")

    # ------------------------------------------------------------------
    # Domain resolution priority:
    # 1. country_name
    # 2. region_name
    # 3. bbox_bounds
    # ------------------------------------------------------------------
    if country_name is not None:
        bbox_bounds = grid_functions.get_country_bbox_bounds(country_name, resolution="110m")
        region_name = None
    elif region_name is not None:
        cname = region_name.lower().replace(" ", "")
        if cname not in grid_functions.DOMAIN_BOUNDS:
            raise ValueError(
                f"Unsupported region_name '{region_name}'. "
                f"Choose from {list(grid_functions.DOMAIN_BOUNDS)}."
            )
        bbox_bounds = grid_functions.DOMAIN_BOUNDS[cname]
        region_name = cname
    else:
        cname = grid_functions.match_domain_from_bbox_bounds(bbox_bounds)
        if cname is not None:
            region_name = cname

    # Snap bbox to the nearest valid global model grid while ensuring
    # the requested domain remains fully contained.
    original_bbox_bounds = list(map(float, bbox_bounds))
    bbox_bounds = grid_functions.snap_bbox_to_grid_containing(bbox_bounds, grid=grid)
    bbox_bounds = [float(x) for x in bbox_bounds]

    # Date handling with rollback to get the latest available forecast
    if fcdate is None:
        fcdate = datetime.datetime.utcnow().strftime("%Y%m%d")

    # Make sure only uppercase model name is used if the provided model contains lowercase.
    model = model.upper() if any(c.islower() for c in model) else model
    origin_id = argument_output.output_originID(model, fcdate)

    # Some forecast dates may not be valid for the given model/origin.
    # In that case, roll back one day at a time until a valid date is found.
    while True:
        try:
            argument_check.check_fcdate(fcdate, origin_id)
            break
        except ValueError:
            old_date = fcdate
            fcdate = (
                datetime.datetime.strptime(fcdate, "%Y%m%d")
                - datetime.timedelta(days=1)
            ).strftime("%Y%m%d")
            if verbose:
                print(f"[INFO] {old_date} not valid, rolling back to {fcdate}...")
   
    # convert start_lt to leadtime if dayofweek
    start_lt = argument_output.convert_dayofweek_to_leadtime(fcdate, start_lt)

    # if leadtime_hour is set, no aggregation!
    if leadtime_hour:
        aggregation_switch = False
    else:
        aggregation_switch = True

    # Get parameters
    leveltype, plevs, webapi_param, ecds_varname, origin_id, leadtime_hour, fc_enslags = (
        argument_output.check_and_output_all_fc_arguments(
            variable,
            model,
            fcdate,
            bbox_bounds,
            data_format,
            grid_for_request,
            plevs,
            leadtime_hour,
            start_lt,
            end_lt,
            period,
            fc_enslags
        )
    )

    # Filename construction
    if filename is None:
        plev_str = ""
        if plevs is not None:
            plevs = [plevs] if isinstance(plevs, int) else plevs
            if len(plevs) == 1:
                plev_str = f"_{plevs[0]}hPa"
            else:
                plev_str = f"_{plevs[0]}-{plevs[-1]}hPa"

        if country_name:
            filename = f"{variable}_{model}_{fcdate}{plev_str}_{country_name.strip().lower().replace(' ', '')}_fc"
        elif region_name:
            filename = f"{variable}_{model}_{fcdate}{plev_str}_{region_name}_fc"
        else:
            north, west, south, east = bbox_bounds
            bounds_str = (
                f"{format_coord(north, lat=True)}_{format_coord(west, lat=False)}_"
                f"{format_coord(south, lat=True)}_{format_coord(east, lat=False)}"
            )
            filename = f"{variable}_{model}_{fcdate}{plev_str}_{bounds_str}_fc"

    # Add extension
    filename_save = filename
    ext = ".nc" if data_format.lower() == "netcdf" else ".grib"
    if not filename_save.endswith(ext):
        filename_save = f"{filename_save}{ext}"

    # Ensure save directory
    if data_save_dir is not None:
        os.makedirs(data_save_dir, exist_ok=True)
        filename_save = os.path.join(data_save_dir, os.path.basename(filename_save))

    print("[DEBUG] Forecast download request")
    print(f"   Model               : {model}")
    print(f"   Target grid         : {grid}")
    print(f"   Variable            : {variable}")
    print(f"   Pressure levels     : {plevs}")
    print(f"   Forecast date       : {fcdate}")

    if country_name is not None:
        print(f"   Requested country   : {country_name}")
    elif region_name is not None:
        print(f"   Requested region    : {region_name}")

    print(f"   Requested bbox      : {original_bbox_bounds}")
    print(f"   Selected bbox       : {bbox_bounds}")
    print(f"   Output file         : {filename_save}")

    # Skip download if file exists
    if os.path.exists(filename_save) and not overwrite:
        print(f"[INFO] File already exists: {filename_save}, skipping download.")
        return filename_save

    # Downstream request code expects path without the .nc extension at this stage
    filename__save = filename_save
    if filename_save.endswith(".nc"):
        filename_save = filename_save[:-3]

    # Print clean info
    if verbose:
        print(f"[INFO] Downloading {variable} (plevs={plevs}) for {country_name or region_name or bbox_bounds}")
        print(f"[INFO] Saving as {filename_save}")

    # ECDS REQUEST FROM S2S DATABASE
    if variable not in ['TC_TRACKS','MJO']:
        try:
            if verbose:
                ecdsAPI_requests.request_forecast(fcdate,origin_id,grid_for_request,variable,bbox_bounds,
                                                  data_format,webapi_param,leadtime_hour,period,
                                                  leveltype,filename_save,plevs,fc_enslags,start_lt,aggregation_switch,cleanup=cleanup)
            else:
                with SuppressOutput():
                    ecdsAPI_requests.request_forecast(fcdate,origin_id,grid_for_request,variable,bbox_bounds,
                                                  data_format,webapi_param,leadtime_hour,period,
                                                  leveltype,filename_save,plevs,fc_enslags,start_lt,aggregation_switch,cleanup=cleanup)
        except Exception:
            print(f"[ERROR] Download failed for {filename_save}")
            raise

    elif variable == 'TC_TRACKS':
        try:
            if verbose:
                download_S2Stc_tracks.download_forecast_TCtracks(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags)
            else:
                with SuppressOutput():
                    download_S2Stc_tracks.download_forecast_TCtracks(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags)
        except Exception:
            print(f"[ERROR] Download failed for {filename_save}")
            raise
    elif variable == 'MJO':
        try:
            if verbose:
                download_S2SMJO_indices.download_forecast_MJO(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags)
            else:
                with SuppressOutput():
                    download_S2SMJO_indices.download_forecast_MJO(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags)
        except Exception:
            print(f"[ERROR] Download failed for {filename_save}")
            raise

    return filename_save
