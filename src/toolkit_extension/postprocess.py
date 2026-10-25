# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

import xarray as xr
import numpy as np
import pandas as pd
import json
import os
from contextlib import contextmanager

import xarray_regrid  # noqa: F401  (registers the .regrid accessor)
import xarray_regrid.methods.conservative as _xrc



#this determines with functions will be publicly visible in the installed package
__all__=[
         "read_netcdf",
         "align_grid",
         "set_verbose",
         "aggregate",
         "hcst_to_effective_lead_time",
         "fcst_to_effective_lead_time",
         "align_obs_to_model_time",
         "align_obs_to_product",
        ]

#######################################################################################

#
# helper functions
#
#######################################################################################

VERBOSE = True

def set_verbose(v):
    """
    Switch the package's log messages (_log) on or off.

    Parameters
    ----------
    v : bool
        True to print log messages, False to silence them (messages logged with
        force=True are always printed). The new setting is printed.
    """
    global VERBOSE
    VERBOSE = v
    print("verbosity is {}".format(v))


def _log(msg, force=False):
    """
    Print `msg` if verbosity is on (see set_verbose), or always if force=True.
    """
    if VERBOSE | force:
        print(msg)

def _check_grid_overlap(obs, hindcast, raise_if_missing=True, threshold=0.9):
    """
    Check spatial overlap between two latitude/longitude grids.

    The overlap is the intersection of the two latitude/longitude bounding boxes, as a
    fraction of the smaller box.

    Parameters
    ----------
    obs : xarray.DataArray or xarray.Dataset
        Observational grid with `latitude` and `longitude` coordinates.
    hindcast : xarray.DataArray or xarray.Dataset
        Hindcast/forecast grid with `latitude` and `longitude` coordinates.
    raise_if_missing : bool, default True
        If True, raise a ValueError when the overlap fraction is below `threshold`;
        if False, only log the overlap fraction.
    threshold : float, optional
        Minimum fraction of overlap relative to the smaller domain
        required to proceed (default: 0.9).

    Returns
    -------
    True, unless an error is raised.

    Raises
    ------
    ValueError
        If the grids do not overlap at all (regardless of raise_if_missing), or if the
        overlap fraction is smaller than the threshold and raise_if_missing is True.
    """

    obs_lat_min = float(obs.latitude.min())
    obs_lat_max = float(obs.latitude.max())
    obs_lon_min = float(obs.longitude.min())
    obs_lon_max = float(obs.longitude.max())

    hc_lat_min = float(hindcast.latitude.min())
    hc_lat_max = float(hindcast.latitude.max())
    hc_lon_min = float(hindcast.longitude.min())
    hc_lon_max = float(hindcast.longitude.max())

    # intersection box
    lat_min = max(obs_lat_min, hc_lat_min)
    lat_max = min(obs_lat_max, hc_lat_max)
    lon_min = max(obs_lon_min, hc_lon_min)
    lon_max = min(obs_lon_max, hc_lon_max)

    if lat_min >= lat_max or lon_min >= lon_max:
        raise ValueError("Observation and hindcast grids do not overlap.")

    overlap_area = (lat_max - lat_min) * (lon_max - lon_min)

    obs_area = (obs_lat_max - obs_lat_min) * (obs_lon_max - obs_lon_min)
    hc_area  = (hc_lat_max - hc_lat_min) * (hc_lon_max - hc_lon_min)

    smaller_area = min(obs_area, hc_area)
    overlap_fraction = overlap_area / smaller_area

    if overlap_fraction < threshold:
        if raise_if_missing:
            raise ValueError(f"Grid overlap too small: {overlap_fraction:.2f} "
                             f"(threshold {threshold:.0%})")
        else:
            _log("Grid overlap : {} ".format(np.round(overlap_fraction,2)))
            return True
    else:
        return True


def _harmonize_coords(ds):
    """
    Harmonize coordinate names to the package convention: latitude, longitude and, for
    observation files, time.

    - latitude : coordinates named lat, Latitude, nav_lat or y (also in other letter cases
      that lower-case to lat, nav_lat or y) are renamed to "latitude".
    - longitude: lon, Longitude, nav_lon or x (likewise) are renamed to "longitude".
    - time     : only in files without a "lead_time" coordinate (i.e. observations, not
      forecasts or hindcasts), coordinates named exactly Time, valid_time or date are
      renamed to "time".

    Parameters
    ----------
    ds : xarray.Dataset or xarray.DataArray
        Dataset with potentially non-standard coordinate names.

    Returns
    -------
    xarray.Dataset or xarray.DataArray
        Dataset with coordinates renamed as above.

    Raises
    ------
    ValueError
        If no latitude or longitude coordinate is found after renaming.
    """


    #detecting type of file

    check_time=False
    #forecast and hindcast files have lead time coordinate
    if "lead_time" not in ds.coords:
        check_time=True

    coord_map = {}

    lat_names = ["lat", "Latitude", "nav_lat", "y"]
    lon_names = ["lon", "Longitude", "nav_lon", "x"]
    time_names = ["Time", "valid_time", "date"]

    for name in ds.coords:
        lname = name.lower()

        if name in lat_names or lname in lat_names:
            coord_map[name] = "latitude"

        elif name in lon_names or lname in lon_names:
            coord_map[name] = "longitude"

        if check_time:
            if name in time_names:
                coord_map[name] = "time"

    if len(coord_map)>0:
        _log("Renaming coordinates: {}".format(coord_map))
    else:
        _log("Nothing to do")

    ds = ds.rename(coord_map)

    if "latitude" not in ds.coords:
        raise ValueError("Latitude coordinate not found.")
    if "longitude" not in ds.coords:
        raise ValueError("Longitude coordinate not found.")

    return ds





#######################################################################################
#
# main functions
#
#######################################################################################


def read_netcdf(file, variable, units=None):
    """
    Read a NetCDF file and return the selected variable, checked to be daily.

    The function opens one NetCDF file with ``xarray.open_dataset``, verifies that the
    requested variable exists, renames coordinates to the package convention (latitude,
    longitude and, for observations, time; see _harmonize_coords), optionally sets the
    units, and checks that the data are daily along their time dimension:
      - observations: "time"
      - hindcast    : "lead_time"
      - forecast    : "valid_time" (calendar dates; its lead_time is a 2-D coordinate)
    The first of time, lead_time, valid_time that is a dimension of the variable is
    checked; the median step must be between 23 and 25 hours.

    Parameters
    ----------
    file : str
        Path to a single NetCDF file.
    variable : str
        Name of the variable to extract from the dataset.
    units : str, optional
        Units of the data, e.g. "mm" for daily precipitation totals. If given, it replaces
        the 'units' attribute of the file (files often come with missing or unusual units).
        The values are not converted: give the units the values are actually in.
        The units are carried through aggregation into the products, so for daily totals
        that are summed (e.g. weekly totals), give the units of an amount ("mm"), not of a
        rate ("mm/day").

    Returns
    -------
    xarray.DataArray
        DataArray containing the requested variable.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    OSError
        If the file cannot be opened.
    KeyError
        If the requested variable is not present in the dataset.
    ValueError
        If no latitude/longitude coordinate is found; if the variable has none of the
        dimensions time, lead_time, valid_time; if that dimension has fewer than 3 steps;
        or if its median step is below 23 hours (sub-daily) or above 25 hours
        (multi-day aggregates).
    """

    _log("\n" + "="*60)
    _log(" START: read_netcdf ".center(60, "="))
    _log("="*60 + "\n")

    _log("reading {} from {}".format(variable,file))


    # Check file existence (for simple paths; glob patterns handled by xarray)
    if isinstance(file, str) and not any(c in file for c in "*?[]"):
        if not os.path.exists(file):
            raise FileNotFoundError(f"data file not found: {file}")

    try:
        ds = xr.open_dataset(file)
        _log("success")
    except Exception as e:
        raise OSError("Unable to open dataset: {}".format(file)) from e

    _log("looking for {}".format(variable))

    if variable not in ds:
        raise KeyError("Variable '{}' not found in dataset.".format(variable))
    else:
        _log("found")

    # Make sure lat,lon and time are there, rename if necessary
    _log("Harmonizing coordinates...")
    ds = _harmonize_coords(ds)

    #and we return a data array rather than a dataset
    da = ds[variable]

    # units given by the user replace those of the file
    if units is not None:
        _log("setting units to '{}' (in file: '{}')".format(units, da.attrs.get("units", "")))
        da = da.assign_attrs(units=units)

    _log("check if data are daily...")
    # forecast "valid_time" (its lead_time is a 2-D coordinate, not a dimension)
    for timevar in ["time", "lead_time", "valid_time"]:
        if timevar in da.dims:
            break
    else:
        raise ValueError(f"no time dimension found in {da.dims}; expected one of "
                         "'time', 'lead_time' or 'valid_time'")

    if da.sizes[timevar] > 2:
        dt = pd.Series(da[timevar].values).diff().dropna()   # works for dates and timedeltas

        # dominant timestep
        median_dt = dt.median()

        if median_dt < pd.Timedelta(hours=23):
            raise ValueError(
                f"Dataset appears to be sub-daily (median {timevar} step {median_dt} < 1 day). "
                "This only processes daily data")
        elif median_dt > pd.Timedelta(hours=25):
            raise ValueError(
                f"Dataset appears to be a multi-day aggregate (median {timevar} step {median_dt} > 1 day). "
                "This only processes daily data")

        _log("daily, {} from {} to {}".format(timevar, str(da[timevar].values[0]),
                                              str(da[timevar].values[-1])))
    else:
        raise ValueError(f"Dataset has only {da.sizes[timevar]} {timevar} step(s); "
                         "expected a daily series.")
    _log("all done")
    _log("\n" + "="*60+"\n")

    return da




LAT_NAMES = ("lat", "latitude")
LON_NAMES = ("lon", "longitude")


def _spatial_names(da):
    """Names of the latitude and longitude dimensions of `da`."""
    lat = next((n for n in LAT_NAMES if n in da.dims), None)
    lon = next((n for n in LON_NAMES if n in da.dims), None)
    if lat is None or lon is None:
        raise ValueError(f"no lat/lon dimensions found in {da.dims}; expected one of "
                         f"{LAT_NAMES} and one of {LON_NAMES}")
    return lat, lon


def align_grid(finegrid, coarsegrid, direction="fine_to_coarse", method="conservative",
               raise_if_missing=True, threshold=0.9, nan_threshold=0.5):
    """
    Put two datasets on the same spatial grid, keeping all other dimensions.

    Parameters
    ----------
    finegrid : DataArray, high resolution, e.g. observations (time, latitude, longitude).
    coarsegrid : DataArray, low resolution, e.g. model data
        (valid_time, effective_lead_time, member, latitude, longitude).
        Spatial dimensions may be named lat/lon or latitude/longitude.
    direction : "fine_to_coarse" (default) or "coarse_to_fine".
        "fine_to_coarse": finegrid is regridded to the coarse grid, conservatively
                          (area-weighted mean of the fine cells in each coarse cell).
        "coarse_to_fine": coarsegrid is interpolated bilinearly to the fine grid; fine cells
                          outside the area spanned by the coarse cell centres get the value
                          of the nearest coarse cell.
    method : "conservative" for fine_to_coarse, "bilinear" for coarse_to_fine.
    raise_if_missing, threshold : passed to _check_grid_overlap (grids that do not overlap
        at all always raise an error).
    nan_threshold : fine_to_coarse only, default 0.5. Largest fraction of NaN (by area) in
        a coarse cell for which a value is still computed from the valid part; above it the
        coarse cell is NaN. Same meaning as xESMF's na_thres.

    Returns
    -------
    finealigned, coarsealigned : both on the chosen grid, with the dimension names of that
        grid. Cells without data in either dataset (at any time / lead / member) are NaN in
        both. Attributes are kept, and the regridding settings are added as attributes
        (regrid_direction, regrid_method, regrid_nan_threshold), so that they reach the
        products and can be repeated for new observations (align_obs_to_product).

    Raises
    ------
    ValueError
        If direction or method is not one of the allowed values, if no lat/lon dimensions
        are found, or if the grids do not overlap enough (see _check_grid_overlap).
    """
    _log(" align_grid ".center(60, "="))

    allowed = {"fine_to_coarse": ["conservative"], "coarse_to_fine": ["bilinear"]}
    if direction not in allowed:
        raise ValueError("direction must be 'fine_to_coarse' or 'coarse_to_fine'")
    if method not in allowed[direction]:
        raise ValueError(f"only {', '.join(allowed[direction])} available for {direction}")

    _check_grid_overlap(finegrid, coarsegrid, raise_if_missing, threshold)

    flat, flon = _spatial_names(finegrid)
    clat, clon = _spatial_names(coarsegrid)

    if direction == "fine_to_coarse":
        _log("Regridding finegrid -> coarsegrid (conservative)...")
        # xarray-regrid matches grids by coordinate name: use the coarse grid's names
        src = finegrid.rename({flat: clat, flon: clon})
        target = xr.Dataset(coords={clat: coarsegrid[clat], clon: coarsegrid[clon]})
        finealigned = src.regrid.conservative(target, latitude_coord=clat, time_dim=None,
                                              skipna=True, nan_threshold=nan_threshold)
        finealigned = finealigned.assign_attrs(finegrid.attrs)
        coarsealigned = coarsegrid
        lat, lon = clat, clon
    else:
        _log("Regridding coarsegrid -> finegrid (bilinear)...")
        src = coarsegrid.rename({clat: flat, clon: flon})
        target = xr.Dataset(coords={flat: finegrid[flat], flon: finegrid[flon]})
        lin = src.regrid.linear(target, time_dim=None)
        nearest = src.reindex({flat: finegrid[flat], flon: finegrid[flon]}, method="nearest")
        coarsealigned = lin.fillna(nearest).assign_attrs(coarsegrid.attrs)
        finealigned = finegrid
        lat, lon = flat, flon

    # common spatial mask: cells with data in both datasets (at any time / lead / member)
    def has_data(da):
        other = [d for d in da.dims if d not in (lat, lon)]
        return da.notnull().any(other) if other else da.notnull()

    mask = has_data(finealigned) & has_data(coarsealigned)
    finealigned = finealigned.where(mask)
    coarsealigned = coarsealigned.where(mask)

    # record the settings (where() keeps the attributes)
    settings = dict(regrid_direction=direction, regrid_method=method,
                    regrid_nan_threshold=float(nan_threshold))
    finealigned = finealigned.assign_attrs(settings)
    coarsealigned = coarsealigned.assign_attrs(settings)

    _log(f"done: {int(mask.sum())} of {mask.size} grid cells with data in both")
    _log("=" * 60)
    return finealigned, coarsealigned





# registry of available aggregation functions
AGG_FUNCTIONS = {}

def register_agg(name):
    """Decorator to register an aggregation function under a name."""
    def wrapper(func):
        AGG_FUNCTIONS[name] = func
        return func
    return wrapper


# --- built-ins ---

@register_agg("mean")
def _mean(da, dim, **kwargs):
    return da.mean(dim=dim)

@register_agg("sum")
def _sum(da, dim, **kwargs):
    return da.sum(dim=dim)

@register_agg("count_zero")
def _count_zero(da, dim, **kwargs):
    return (da == 0).sum(dim=dim)

@register_agg("count_above")
def _count_above(da, dim, threshold=0, **kwargs):
    return (da > threshold).sum(dim=dim)

@register_agg("max")
def _max(da, dim, **kwargs):
    return da.max(dim=dim)

@register_agg("max_dry_spell")
def _max_dry_spell(da, dim, threshold=1.0, **kwargs):
    """Longest run of consecutive days with values below `threshold` (dry days) within the
    window. Missing days interrupt a run; NaN if all days are missing."""
    dry = da < threshold                       # NaN compares as False: not dry
    run = xr.zeros_like(dry.isel({dim: 0}, drop=True), dtype=int)
    longest = run
    for i in range(da.sizes[dim]):
        day_dry = dry.isel({dim: i}, drop=True)
        run = (run + 1) * day_dry              # one more day, or back to 0 if not dry
        longest = np.maximum(longest, run)
    return longest.where(da.notnull().any(dim))

# aggregations whose result is a number of days, not a value in the units of the data
COUNT_METHODS = ("count_zero", "count_above", "max_dry_spell")




def aggregate(
    data,
    window_size,
    method="mean",
    dim="effective_lead_time",
    **kwargs):
    """
    Aggregate daily data over `dim` in consecutive windows of `window_size` days.

    The lead times are converted to timedeltas (integers are taken as days) and resampled
    into windows of `window_size` days with label="left", closed="left": each output lead
    time is the start of its window, e.g. with window_size=7, lead 7 days covers the daily
    leads 7 to 13. If the number of steps along `dim` is not a multiple of window_size,
    the incomplete trailing window is dropped (with a log message). Grid cells that
    are NaN at all leads (and members / valid times) in the input are NaN in the output.

    Parameters
    ----------
    data : xr.DataArray or xr.Dataset
        Daily data, e.g. (valid_time, effective_lead_time, [member,] latitude, longitude).
    window_size : int
        Window length in days, between 5 and 50.
    method : str, default "mean"
        Name of a function registered in AGG_FUNCTIONS: "mean", "sum", "max",
        "count_zero" (number of days equal to 0), "count_above" (number of days above
        `threshold`, passed as a keyword argument; default 0), "max_dry_spell" (longest
        run of consecutive days below `threshold` within the window; default 1.0, e.g.
        1 mm for precipitation; a spell that crosses into the next window is split), or
        one added with register_agg.
    dim : str, default "effective_lead_time"
        Dimension being aggregated.
    **kwargs : passed through to the aggregation function (e.g. threshold=1.0).

    Returns
    -------
    Same type as data, with `dim` at the window starts. The attributes of the input are
    kept and the settings are added (aggregation_method, aggregation_window and, if any,
    aggregation_kwargs as JSON text), so that they reach the products and can be repeated
    for new observations (align_obs_to_product). For counting methods (COUNT_METHODS) the
    units are set to "days".

    Raises
    ------
    ValueError
        If method is not registered or window_size is not between 5 and 50.
    """

    _log(f"Aggregating lead_time to {window_size}-D blocks using ({method})...")

    #arbitrary choice
    agg_windows=list(range(5,51))

    #checking if method is available
    if method not in AGG_FUNCTIONS:
        raise ValueError(
            f"Unknown aggregation method '{method}'. "
            f"Available: {list(AGG_FUNCTIONS)}"
        )

    #picking up requested method
    func = AGG_FUNCTIONS[method]

    # validating aggregation window
    if window_size not in agg_windows:
        raise ValueError(f"window_size must be between 5 and 50, got {window_size}")

    #agg_window is just window in the xarray format
    agg_window=f"{window_size}D"

    # this converts integer index to timedelta - and this is necessary for the aggregator function
    data = data.assign_coords({dim: pd.to_timedelta(data[dim], "D")})

    # truncate to a multiple of window_days before resampling
    n_total = data.sizes[dim]
    n_full = (n_total // window_size) * window_size

    if n_full < n_total:
        dropped = n_total - n_full
        _log(f"aggregate: dropping last {dropped} step(s) of '{dim}' ({n_total} not divisible "
             f"by window_size={window_size}); keeping {n_full} steps.")
        data = data.isel({dim: slice(0, n_full)})

    resampler = data.resample(
        **{dim: agg_window},
        label="left",
        closed="left"
    )

    #aggregating done here!
    aggdata = resampler.map(lambda da: func(da, dim=dim, **kwargs))


    #finding nans
    dims=[dim]
    for v in ["member","target_time","valid_time","init_date","init_lag"]:
        if v in data.dims:
            dims.append(v)

    spatial_nan_mask = np.isnan(data).all(dims)

    #removing all nans
    aggdata = aggdata.where(~spatial_nan_mask)

    # attributes of the input plus the aggregation settings
    attrs = dict(data.attrs)
    attrs.update(aggregation_method=method, aggregation_window=int(window_size))
    if kwargs:
        attrs["aggregation_kwargs"] = json.dumps(kwargs)
    if method in COUNT_METHODS:
        attrs["units"] = "days"
    aggdata = aggdata.assign_attrs(attrs)

    return aggdata








def _find_cycles(init_dates, max_gap=None):
    """
    Group sorted init dates into cycles of staggered initialisations.

    A new cycle starts wherever the gap to the previous init date exceeds `max_gap`.
    max_gap : timedelta-like (e.g. "10D"), or None to infer it: with all gaps equal there is
              a single cycle (e.g. a forecast); otherwise the threshold is the geometric mean
              of the smallest and largest gap (2 days and ~1 year -> ~27 days).
    Returns a 2-D array (n_cycles, n_inits_per_cycle) of dates; all cycles must have the
    same number of inits.
    """
    dates = np.sort(np.asarray(init_dates))
    if len(dates) == 1:
        return dates[None, :]
    gaps = np.diff(dates)
    if max_gap is None:
        g = gaps / np.timedelta64(1, "s")
        if g.min() == g.max():
            return dates[None, :]                               # one cycle
        thresh = np.timedelta64(int(np.sqrt(g.min() * g.max())), "s")
    else:
        thresh = np.timedelta64(pd.Timedelta(max_gap).value, "ns")
    starts = np.r_[0, np.flatnonzero(gaps > thresh) + 1]
    groups = np.split(dates, starts[1:])
    sizes = sorted({len(g) for g in groups})
    if len(sizes) > 1:
        raise ValueError(f"cycles have different numbers of inits {sizes} with a gap threshold "
                         f"of {thresh.astype('timedelta64[h]')}; check the init dates or pass max_gap")
    return np.stack(groups)


def _check_daily_leads(leads, name):
    """Raise unless lead times are whole days, one day apart."""
    leads = np.asarray(leads).astype("timedelta64[ns]")
    steps = np.diff(leads)
    if (leads % ONE_DAY).any() or (len(steps) and (steps != ONE_DAY).any()):
        found = sorted(set((steps / ONE_DAY).round(3).tolist())) if len(steps) else []
        found = [int(f) if float(f).is_integer() else f for f in found]
        raise ValueError(f"'{name}' must be daily (steps of 1 day); found steps of {found} days. "
                         "Restructure daily data first and aggregate afterwards.")


def hcst_to_effective_lead_time(ds, n_inits_per_cycle=None, *, anchor=None, trim=False, max_gap=None,
                      init_dim="hc_init_date", lead_dim="lead_time", member_dim="member",
                      drop_coords=("valid_time",)):
    """
    Combine staggered hindcast initialisations into cycles with a common lead-time axis.

    Parameters
    ----------
    ds : Dataset or DataArray with dims (init_dim, lead_dim, member_dim, ...), e.g.
         (hc_init_date, lead_time, member, latitude, longitude).
         lead_dim must be daily (timedelta, steps of 1 day); otherwise ValueError.
         Aggregate only after restructuring.
    n_inits_per_cycle : number of initialisations per cycle. If None (default) the cycles are
         inferred from gaps between init dates (see _find_cycles); if given, it is checked
         against the inferred grouping.
    anchor : "earliest" or "latest". Date that labels the cycle (valid_time) and from which
         effective_lead_time is counted. Default: "earliest", or "latest" if trim=True.
         "earliest": later inits start later -> NaN at their first leads.
         "latest"  : earlier inits are cut before the anchor and end earlier ->
                     NaN at their last leads.
         The effective_lead_time axis has the values of lead_dim (minus the trailing
         leads removed by trim); it starts at 0 only if lead_dim starts at 0.
    trim : if True, keep only effective lead times covered by every initialisation (no NaN
         from the staggering). Requires anchor="latest" (the default when trim=True); the
         trailing leads where earlier inits have already ended are removed.
    max_gap : passed to _find_cycles; only needed if the automatic grouping fails.
    init_dim, lead_dim, member_dim : names of the init, lead and member dimensions of ds.
    drop_coords : coordinates along (init_dim, lead_dim) to drop; they would be wrong after
         shifting. 'target_time' (= valid_time + effective_lead_time) is added instead.

    Returns
    -------
    Same type as ds, dims (valid_time, effective_lead_time, member, ...), with
      valid_time = the anchor init date of each cycle, effective_lead_time timedelta64[ns],
      member = members of all inits of the cycle (n_inits x n_members), and coordinates
      init_date(valid_time, member) : actual initialisation date of each member
      init_lag(valid_time, member)  : init_date - valid_time (timedelta)
      target_time(valid_time, effective_lead_time)
    plus any other coordinates along (init_dim, member_dim), e.g. 'lag', as (valid_time, member).
    Data variables are float with NaN where an initialisation has no data; their on-disk
    packing (integer dtype, scale_factor, _FillValue) is cleared so that NaN survives saving.

    Raises
    ------
    ValueError
        If anchor is invalid or trim=True with anchor="earliest"; if lead_dim is not daily;
        if the cycles found do not match n_inits_per_cycle or have different numbers of
        inits; if inits within a cycle are not whole days apart.
    """
    if anchor is None:
        anchor = "latest" if trim else "earliest"
    if anchor not in ("earliest", "latest"):
        raise ValueError("anchor must be 'earliest' or 'latest'")
    if trim and anchor != "latest":
        raise ValueError("trim=True requires anchor='latest'")
    _check_daily_leads(ds[lead_dim].values, lead_dim)
    ds = ds.drop_vars([c for c in drop_coords if c in ds.coords]).sortby(init_dim)

    cycles = _find_cycles(ds[init_dim].values, max_gap)
    if n_inits_per_cycle is not None and cycles.shape[1] != n_inits_per_cycle:
        raise ValueError(f"found {cycles.shape[1]} inits per cycle in the data, "
                         f"but n_inits_per_cycle={n_inits_per_cycle}")
    if ((cycles - cycles[:, :1]) % ONE_DAY).any():
        raise ValueError("inits within a cycle must be whole days apart")

    leads = ds[lead_dim].values
    n_mem = ds.sizes[member_dim]

    out_cycles = []
    for dates in cycles:
        a = dates[0] if anchor == "earliest" else dates[-1]
        parts = []
        for k, d in enumerate(dates):
            offset = d - a                                   # >= 0 (earliest) or <= 0 (latest)
            x = ds.sel({init_dim: d}, drop=True)
            # value at original lead L belongs to effective lead L + offset
            x = x.assign_coords({lead_dim: x[lead_dim] + offset}).reindex({lead_dim: leads})
            x = x.assign_coords({member_dim: np.arange(k * n_mem, (k + 1) * n_mem,
                                                       dtype=ds[member_dim].dtype),
                                 "init_date": (member_dim, np.full(n_mem, d)),
                                 "init_lag": (member_dim, np.full(n_mem, offset))})
            parts.append(x)
        cyc = xr.concat(parts, dim=member_dim, coords="minimal", compat="equals", join="exact")
        out_cycles.append(cyc.expand_dims(valid_time=[a]))

    member_coords = [c for c in out_cycles[0].coords
                     if member_dim in out_cycles[0][c].dims and c != member_dim]
    out = xr.concat(out_cycles, dim="valid_time", coords=member_coords,
                    compat="equals", join="exact")
    out = out.rename({lead_dim: "effective_lead_time"})
    out = out.assign_coords(effective_lead_time=out["effective_lead_time"].astype("timedelta64[ns]"))

    if trim:
        # leads covered by every init in every cycle, from the offsets (not from the data,
        # which may contain genuine NaN, e.g. over the ocean)
        lags = out["init_lag"].values                        # (valid_time, member), <= 0
        last = leads[-1] + lags.min()                        # earliest init ends first
        out = out.sel(effective_lead_time=slice(leads[0], last))

    out = out.assign_coords(target_time=out["valid_time"] + out["effective_lead_time"])

    _clear_packing(out)
    return out.transpose("valid_time", "effective_lead_time", member_dim, ...)






def fcst_to_effective_lead_time(fc, *, leads=None, time_dim="valid_time", member_dim="member",
                                init_coord="fc_init_date", drop_coords=("lead_time",)):
    """
    Bring a forecast that is already on a calendar-date axis into the same format as
    hcst_to_effective_lead_time(hindcast, trim=True).

    Input  : (valid_time, member, latitude, longitude), daily, where valid_time are target
             dates and each member has its initialisation date in coordinate `init_coord`
             (member,).
    Output : (valid_time: 1, effective_lead_time, member, latitude, longitude)
             valid_time          = first date (from the latest initialisation on) at which
                                   every member has data
             effective_lead_time = target date - valid_time, starting at 0. This matches
                                   the hindcast only if the hindcast's lead_time starts at 0
                                   (hcst_to_effective_lead_time keeps the lead_time values).
             coords init_date, init_lag (= init_date - valid_time) (valid_time, member),
             target_time (valid_time, effective_lead_time), other member coords (e.g. 'lag')
             as (valid_time, member).

    Parameters
    ----------
    fc : Dataset or DataArray, as above.
    leads : effective lead times to keep, e.g. hcst.effective_lead_time, so that forecast and
            hindcast have identical lead axes; an error is raised if a requested lead is
            missing or if any member has no data at it. If None, all leads up to the last one
            at which every member has data are kept (like trim=True).
    time_dim : name of the calendar-date dimension of fc (default "valid_time").
    member_dim : name of the member dimension (default "member").
    init_coord : name of the coordinate with each member's initialisation date (default
            "fc_init_date"); it is renamed to init_date.
    drop_coords : coordinates to drop, e.g. the forecast's own lead_time.

    Raises
    ------
    ValueError
        If time_dim is not daily; if no date has data from every member; if `leads` are
        missing or have members without data.
    """
    t = fc[time_dim].values
    _check_daily_leads(t - t[0], time_dim)
    fc = fc.drop_vars([c for c in drop_coords if c in fc.coords])
    anchor = fc[init_coord].values.max()

    def gaps(o):
        """True at dates where any member has no data anywhere in the domain."""
        first = o if isinstance(o, xr.DataArray) else o[list(o.data_vars)[0]]
        space = [d for d in first.dims if d not in (time_dim, member_dim)]
        g = first.isnull().all(space) if space else first.isnull()
        return g.any(member_dim).values

    # dates from the latest init on; the first date with data from every member is lead 0
    out = fc.isel({time_dim: np.flatnonzero(fc[time_dim].values >= anchor)})
    ok = np.flatnonzero(~gaps(out))
    if ok.size == 0:
        raise ValueError("no date at which every member has data")
    out = out.isel({time_dim: slice(ok[0], None)})
    dates = out[time_dim].values
    out = (out.assign_coords(effective_lead_time=(time_dim, (dates - dates[0]).astype("timedelta64[ns]")),
                             target_time=(time_dim, dates))
              .swap_dims({time_dim: "effective_lead_time"})
              .drop_vars(time_dim))

    def gap_leads(o):
        return gaps(o.rename({"effective_lead_time": time_dim}))

    if leads is None:
        bad = np.flatnonzero(gap_leads(out))
        if bad.size:                                       # trailing gaps: earlier inits ended
            out = out.isel(effective_lead_time=slice(0, bad[0]))
    else:
        leads = np.asarray(leads).astype("timedelta64[ns]")
        missing = np.setdiff1d(leads, out["effective_lead_time"].values)
        if missing.size:
            raise ValueError(f"forecast does not cover leads {missing.astype('timedelta64[D]')}")
        out = out.sel(effective_lead_time=leads)
        bad = out["effective_lead_time"].values[gap_leads(out)]
        if bad.size:
            raise ValueError(f"some members have no data at leads {bad.astype('timedelta64[D]')}")

    vt = dates[0]                                          # effective lead 0
    out = out.drop_vars("target_time")
    out = out.rename({init_coord: "init_date"}) if init_coord != "init_date" else out
    out = out.assign_coords(init_lag=out["init_date"] - vt)
    out = out.expand_dims(valid_time=[vt])
    member_coords = [c for c in out.coords if out[c].dims == (member_dim,) and c != member_dim]
    out = out.assign_coords({c: out[c].expand_dims(valid_time=[vt]) for c in member_coords})
    out = out.assign_coords(target_time=out["valid_time"] + out["effective_lead_time"])
    _clear_packing(out)
    return out.transpose("valid_time", "effective_lead_time", member_dim, ...)




def _clear_packing(obj):
    """NaN must survive writing: drop integer packing inherited from the source file."""
    variables = [obj] if isinstance(obj, xr.DataArray) else [obj[v] for v in obj.data_vars]
    for var in variables:
        for key in ("dtype", "scale_factor", "add_offset", "_FillValue", "missing_value"):
            var.encoding.pop(key, None)


ONE_DAY = np.timedelta64(1, "D")




def align_obs_to_model_time(obs, model, *, add_members=False, raise_if_missing=True, obs_time="time",
              init_dim="valid_time", lead_dim="effective_lead_time", member_dim="member"):
    """
    Put daily observations on the (valid_time, effective_lead_time) structure of `model`.

    The obs value at target date valid_time + effective_lead_time is taken for every
    (valid_time, lead) of the model. Dates of both obs and model are normalised to
    midnight first.

    Parameters
    ----------
    obs : DataArray or Dataset (time, latitude, longitude), daily, no gaps.
    model : hindcast or forecast in effective-init format, with daily lead times.
    add_members : replicate obs along the model's member dim, NaN where a model member has
          no data anywhere in the domain (e.g. staggering gaps in untrimmed data).
    raise_if_missing : if True, raise if obs do not cover all model target dates; if False,
          drop the valid_times (cycles) that are not fully covered, from obs and model.
    obs_time : name of the time dimension of obs (default "time").
    init_dim, lead_dim, member_dim : names of the model's valid_time, lead and member
          dimensions (defaults "valid_time", "effective_lead_time", "member").

    Returns
    -------
    obs_aligned : (valid_time, effective_lead_time, [member,] latitude, longitude), with
          target_time if the model has that coordinate.
    model_aligned : model with valid_time normalised to midnight, restricted to fully
          covered valid_times if raise_if_missing=False.

    Raises
    ------
    ValueError if obs are not daily without gaps, if the model lead times are not daily
    (e.g. already aggregated), or (raise_if_missing=True) if obs do not cover all model
    target dates.
    """
    _log(" align_obs ".center(60, "="))

    obs = obs.assign_coords({obs_time: obs.indexes[obs_time].normalize()})
    model = model.assign_coords({init_dim: model.indexes[init_dim].normalize()})

    if len(obs[obs_time]) > 1 and (np.diff(obs[obs_time].values) != ONE_DAY).any():
        raise ValueError("obs must be daily without gaps")
    leads = model[lead_dim].values.astype("timedelta64[ns]")
    if len(leads) > 1 and (np.diff(leads) != ONE_DAY).any():
        raise ValueError(f"model lead times must be daily; found steps of "
                         f"{sorted(set(np.diff(leads).astype('timedelta64[D]').astype(int).tolist()))} "
                         "days. Align daily data first, then aggregate model and obs together.")

    # target dates (valid_time, lead)
    vt = model[init_dim].values
    dates = vt[:, None] + leads[None, :]

    # coverage
    have = obs[obs_time].values
    missing_per_vt = (~np.isin(dates, have)).any(axis=1)
    _log(f"obs: {str(have[0])[:10]} to {str(have[-1])[:10]} | "
         f"model targets: {str(dates.min())[:10]} to {str(dates.max())[:10]}")
    if missing_per_vt.any():
        if raise_if_missing:
            raise ValueError(f"obs do not cover {missing_per_vt.sum()} of {len(vt)} valid_times")
        _log(f"WARNING: dropping {missing_per_vt.sum()} valid_times not covered by obs")
        keep = ~missing_per_vt
        model, dates, vt = model.isel({init_dim: keep}), dates[keep], vt[keep]
    else:
        _log("obs cover all model target dates")

    # select obs for every (valid_time, lead) at once
    out = obs.sel({obs_time: xr.DataArray(dates, dims=[init_dim, lead_dim])}).drop_vars(obs_time)
    out = out.assign_coords({init_dim: vt, lead_dim: model[lead_dim].values})
    if "target_time" in model.coords:
        out = out.assign_coords(target_time=model["target_time"])

    if add_members:
        first = model if isinstance(model, xr.DataArray) else model[list(model.data_vars)[0]]
        space = [d for d in first.dims if d not in (init_dim, lead_dim, member_dim)]
        has_data = first.notnull().any(space)                 # (valid_time, lead, member)
        member_coords = {c: model[c] for c in model.coords
                         if member_dim in model[c].dims and c != member_dim}
        out = out.expand_dims({member_dim: model[member_dim].values}).where(has_data)
        out = out.assign_coords(member_coords)

    order = [init_dim, lead_dim] + ([member_dim] if add_members else [])
    out = out.transpose(*order, ...)
    _log("=" * 60)
    return out, model




def align_obs_to_product(obs, product, *, raise_if_missing=True):
    """
    Put raw daily observations on the grid and (aggregated) lead times of a saved forecast
    product, repeating the steps of the hindcast pipeline with the settings stored in the
    product, so that the result can go straight into verify_deterministic /
    verify_probabilistic.

    Steps (settings from the product attributes, written there by align_grid and aggregate):
      1. grid  - regrid_direction "fine_to_coarse": obs regridded conservatively to the
                 product grid with align_grid (regrid_method, regrid_nan_threshold);
                 "coarse_to_fine": the obs kept their own grid in the pipeline (the model
                 was interpolated to it), so the product's grid points are taken from the
                 obs, which must be on the same grid.
                 Grid cells without forecast data are set to NaN.
      2. dates - the daily obs of every day covered by the forecast: valid_time + daily
                 leads from the first lead to the end of the last aggregation window
                 (with align_obs_to_model_time).
      3. aggregation - aggregate() with aggregation_method, aggregation_window and
                 aggregation_kwargs; lead times are the starts of the windows, as in the
                 product.

    Parameters
    ----------
    obs : DataArray (time, latitude, longitude), daily, as returned by read_netcdf
          (with units=... if the file's units need replacing).
    product : product Dataset (make_deterministic_product / make_probabilistic_product, or
          the saved file).
    raise_if_missing : passed to align_obs_to_model_time and, for "fine_to_coarse", to
          align_grid; with True (default) an error is raised if the obs do not cover the
          forecast dates or the grid.

    Returns
    -------
    DataArray (valid_time, effective_lead_time, latitude, longitude), aligned with the
    product's forecast.

    Raises
    ------
    ValueError
        If the product lacks the regrid_* / aggregation_* attributes; if (coarse_to_fine)
        the obs are on a different grid; if the obs do not cover the forecast dates; if the
        aggregated leads differ from the product's; if the obs units differ from the
        product's.
    """
    _log(" align_obs_to_product ".center(60, "="))

    a = product.attrs
    needed = ["regrid_direction", "regrid_method", "regrid_nan_threshold",
              "aggregation_method", "aggregation_window"]
    missing = [k for k in needed if k not in a]
    if missing:
        raise ValueError(f"the product has no {missing} attributes: it was made from obs that "
                         "did not go through the current align_grid / aggregate; regenerate it")

    fc = product["fcst_mean"] if "fcst_mean" in product else product["fcst_prob"]
    other = [d for d in fc.dims if d not in ("latitude", "longitude")]
    has_data = fc.notnull().any(other)                    # (latitude, longitude)

    # 1. grid
    if a["regrid_direction"] == "fine_to_coarse":
        template = xr.ones_like(has_data, dtype=float).where(has_data)
        obs_grid, _ = align_grid(obs, template, direction="fine_to_coarse",
                                 method=a["regrid_method"],
                                 nan_threshold=float(a["regrid_nan_threshold"]),
                                 raise_if_missing=raise_if_missing)
    else:
        try:
            obs_grid = obs.sel(latitude=fc["latitude"].values, longitude=fc["longitude"].values)
        except KeyError:
            raise ValueError("the product was made with regrid_direction='coarse_to_fine', i.e. on "
                             "the grid of the obs, but these obs are on a different grid") from None
    obs_grid = obs_grid.drop_vars([c for c in obs_grid.coords if c not in obs_grid.dims])
    obs_grid = obs_grid.where(has_data)

    # 2. daily dates covered by the aggregated leads
    window = int(a["aggregation_window"])
    leads = product["effective_lead_time"].values.astype("timedelta64[ns]")
    daily = np.arange(leads[0], leads[-1] + window * ONE_DAY, ONE_DAY).astype("timedelta64[ns]")
    model_like = xr.DataArray(np.zeros((product.sizes["valid_time"], len(daily))),
                              dims=["valid_time", "effective_lead_time"],
                              coords=dict(valid_time=product["valid_time"].values,
                                          effective_lead_time=daily))
    obs_daily, _ = align_obs_to_model_time(obs_grid, model_like,
                                           raise_if_missing=raise_if_missing)

    # 3. aggregation, as for the hindcast
    kwargs = json.loads(a.get("aggregation_kwargs", "{}"))
    out = aggregate(obs_daily, window, a["aggregation_method"], **kwargs)
    found = out["effective_lead_time"].values.astype("timedelta64[ns]")
    if not np.array_equal(found, leads):
        raise ValueError(f"aggregated obs leads {found.astype('timedelta64[D]').tolist()} differ "
                         f"from the product's {leads.astype('timedelta64[D]').tolist()}; the "
                         "product's aggregation attributes do not match its lead times")
    out = out.assign_coords(effective_lead_time=leads).where(has_data)

    u_obs, u_prod = out.attrs.get("units"), a.get("units")
    if u_obs and u_prod and u_obs != u_prod:
        raise ValueError(f"obs units '{u_obs}' differ from the product's '{u_prod}'; set the "
                         "units when reading the obs (read_netcdf(..., units=...))")

    _log("=" * 60)
    return out.transpose("valid_time", "effective_lead_time", "latitude", "longitude")






@contextmanager
def _dense_weights():
    """Make xarray-regrid use dense weights even if `sparse` is installed.
    Relies on xarray-regrid's module attribute `sparse` (checked in v0.4.2)."""
    saved = getattr(_xrc, "sparse", None)
    _xrc.sparse = None
    try:
        yield
    finally:
        _xrc.sparse = saved


def _crop_to(src, slat, slon, ref, rlat, rlon):
    """Crop src to the extent of ref plus one ref grid cell (either latitude order)."""
    def step(c):
        v = np.asarray(c.values, dtype=float)
        return float(np.abs(np.diff(v)).max()) if v.size > 1 else 0.0
    pad_lat, pad_lon = step(ref[rlat]), step(ref[rlon])
    lo, hi = float(ref[rlat].min()) - pad_lat, float(ref[rlat].max()) + pad_lat
    lat_sl = slice(hi, lo) if src[slat].values[0] > src[slat].values[-1] else slice(lo, hi)
    lon_sl = slice(float(ref[rlon].min()) - pad_lon, float(ref[rlon].max()) + pad_lon)
    out = src.sel({slat: lat_sl, slon: lon_sl})
    if out.sizes[slat] == 0 or out.sizes[slon] == 0:
        raise ValueError("grids do not overlap after cropping; check longitude conventions "
                         "(0..360 vs -180..180) and latitude ranges")
    return out


def align_grid_sparse(finegrid, coarsegrid, direction="fine_to_coarse", method="conservative",
               raise_if_missing=True, threshold=0.9, nan_threshold=0.5):
    """
    Put two datasets on the same spatial grid, keeping all other dimensions; a variant of
    align_grid for environments where the package `sparse` is installed.

    If `sparse` is installed (it is a dependency of xESMF, among others), xarray-regrid
    stores its conservative weights as sparse arrays. For dense data this is ~300x slower,
    uses far more memory, and returns sparse-backed results. This function therefore runs
    xarray-regrid with dense weights; align_grid does not, and is meant for environments
    without `sparse`.

    Differences from align_grid:
      - the source grid is first cropped to the extent of the target grid plus one target
        grid cell (_crop_to), which makes regridding faster;
      - conservative regridding uses dense weights (_dense_weights).
    As in align_grid, the regridding settings are added as attributes, so products made
    from its output can be used with align_obs_to_product (which itself uses align_grid).

    Parameters
    ----------
    finegrid : DataArray, high resolution, e.g. observations (time, latitude, longitude).
    coarsegrid : DataArray, low resolution, e.g. model data
        (valid_time, effective_lead_time, member, latitude, longitude).
        Spatial dimensions may be named lat/lon or latitude/longitude.
    direction : "fine_to_coarse" (default) or "coarse_to_fine".
        "fine_to_coarse": finegrid is regridded to the coarse grid, conservatively
                          (area-weighted mean of the fine cells in each coarse cell).
        "coarse_to_fine": coarsegrid is interpolated bilinearly to the fine grid; fine cells
                          outside the area spanned by the coarse cell centres get the value
                          of the nearest coarse cell.
    method : "conservative" for fine_to_coarse, "bilinear" for coarse_to_fine.
    raise_if_missing, threshold : passed to _check_grid_overlap (grids that do not overlap
        at all always raise an error).
    nan_threshold : fine_to_coarse only, default 0.5. Largest fraction of NaN (by area) in
        a coarse cell for which a value is still computed from the valid part; above it the
        coarse cell is NaN. Same meaning as xESMF's na_thres.

    Returns
    -------
    finealigned, coarsealigned : both on the chosen grid, with the dimension names of that
        grid. Cells without data in either dataset are NaN in both. Attributes are kept,
        and regrid_direction, regrid_method and regrid_nan_threshold are added.

    Raises
    ------
    ValueError
        If direction or method is not one of the allowed values, if no lat/lon dimensions
        are found, or if the grids do not overlap (enough).
    """
    _log(" align_grid ".center(60, "="))

    allowed = {"fine_to_coarse": ["conservative"], "coarse_to_fine": ["bilinear"]}
    if direction not in allowed:
        raise ValueError("direction must be 'fine_to_coarse' or 'coarse_to_fine'")
    if method not in allowed[direction]:
        raise ValueError(f"only {', '.join(allowed[direction])} available for {direction}")

    _check_grid_overlap(finegrid, coarsegrid, raise_if_missing, threshold)

    flat, flon = _spatial_names(finegrid)
    clat, clon = _spatial_names(coarsegrid)

    if direction == "fine_to_coarse":
        _log("Regridding finegrid -> coarsegrid (conservative)...")
        # crop first: xarray-regrid's cost grows with the size of the source grid
        # xarray-regrid matches grids by coordinate name: use the coarse grid's names
        src = _crop_to(finegrid, flat, flon, coarsegrid, clat, clon).rename({flat: clat, flon: clon})
        target = xr.Dataset(coords={clat: coarsegrid[clat], clon: coarsegrid[clon]})
        with _dense_weights():
            finealigned = src.regrid.conservative(target, latitude_coord=clat, time_dim=None,
                                                  skipna=True, nan_threshold=nan_threshold)
        finealigned = finealigned.assign_attrs(finegrid.attrs)
        coarsealigned = coarsegrid
        lat, lon = clat, clon
    else:
        _log("Regridding coarsegrid -> finegrid (bilinear)...")
        src = _crop_to(coarsegrid, clat, clon, finegrid, flat, flon).rename({clat: flat, clon: flon})
        target = xr.Dataset(coords={flat: finegrid[flat], flon: finegrid[flon]})
        lin = src.regrid.linear(target, time_dim=None)
        nearest = src.reindex({flat: finegrid[flat], flon: finegrid[flon]}, method="nearest")
        coarsealigned = lin.fillna(nearest).assign_attrs(coarsegrid.attrs)
        finealigned = finegrid
        lat, lon = flat, flon

    # common spatial mask: cells with data in both datasets (at any time / lead / member)
    def has_data(da):
        other = [d for d in da.dims if d not in (lat, lon)]
        return da.notnull().any(other) if other else da.notnull()

    mask = has_data(finealigned) & has_data(coarsealigned)
    finealigned = finealigned.where(mask)
    coarsealigned = coarsealigned.where(mask)

    # record the settings (where() keeps the attributes)
    settings = dict(regrid_direction=direction, regrid_method=method,
                    regrid_nan_threshold=float(nan_threshold))
    finealigned = finealigned.assign_attrs(settings)
    coarsealigned = coarsealigned.assign_attrs(settings)

    _log(f"done: {int(mask.sum())} of {mask.size} grid cells with data in both")
    _log("=" * 60)
    return finealigned, coarsealigned
