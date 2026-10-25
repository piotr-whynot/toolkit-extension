"""
Forecast products plus the hindcast skill of the same products, one file per forecast.

::

    prob = make_probabilistic_product(obs, hcst, fcst, quantiles=[1/3, 2/3], path="prob.nc")
    rain = make_probabilistic_product(obs, hcst, fcst, thresholds=[1.0], path="rain.nc")
    det  = make_deterministic_product(obs, hcst, fcst, path="det.nc")

Dimensions: valid_time, effective_lead_time ('lead' below), member, latitude and longitude
('lat', 'lon' below).

Inputs (all prepared the same way: restructured, aligned, regridded, aggregated,
bias-corrected; the hindcast cross-validated)::

    obs  : (valid_time, lead, lat, lon)            hindcast years
    hcst : (valid_time, lead, member, lat, lon)
    fcst : (valid_time, lead, member, lat, lon)    usually one valid_time; same leads and grid
                                                   (extra forecast leads are dropped)

Hindcast variables use the dimension 'hcst_valid_time' (the hindcast years), so that they can
sit in the same file as the forecast ('valid_time'); otherwise xarray would merge the two
time axes and fill the gaps with NaN.

Probabilistic product::

    fcst_prob       (valid_time, lead, category, lat, lon)       forecast probabilities
    fcst_n_members   (valid_time, lead, lat, lon)
    hcst_prob        (hcst_valid_time, lead, category, lat, lon)  hindcast probabilities
    hcst_n_members   (hcst_valid_time, lead, lat, lon)
    obs_cat          (hcst_valid_time, lead, category, lat, lon)  observed category (one-hot)
    obs              (hcst_valid_time, lead, lat, lon)            hindcast-period observations
    threshold_obs    (edge, lead, lat, lon)                       to categorise future obs
    threshold_model  (edge, lead, lat, lon)                       only if edges_from="model"
    clim_prob        (category,)                                  reference forecast
                     (lead, category, lat, lon) for absolute thresholds: observed frequency
    skill_<score>    (lead, [category,] lat, lon)                 hindcast skill

Deterministic product (ensemble mean)::

    fcst_mean       (valid_time, lead, lat, lon)
    fcst_anomaly     (valid_time, lead, lat, lon)       fcst_mean - climatology
    fcst_anomaly_pct (valid_time, lead, lat, lon)       100 * fcst_anomaly / climatology
    hcst_mean        (hcst_valid_time, lead, lat, lon)
    hcst_anomaly     (hcst_valid_time, lead, lat, lon)
    hcst_anomaly_pct (hcst_valid_time, lead, lat, lon)
    obs              (hcst_valid_time, lead, lat, lon)  hindcast-period observations
    clim_obs         (lead, lat, lon)                   to compute anomalies of future obs
    clim_model       (lead, lat, lon)                   only if clim_from="model"
    skill_<score>    (lead, lat, lon)                   hindcast skill (of the ensemble mean)

Percent anomalies are meant for positive variables such as precipitation; they are NaN
where the climatology is 0 or negative, and become very large where it is close to 0
(e.g. dry seasons), so mask or cap them for display.

Both carry attributes: product, units, hindcast period, settings, your `attrs`, the
regridding and aggregation settings carried by the obs (regrid_*, aggregation_*, written by
postprocess.align_grid and postprocess.aggregate), which postprocess.align_obs_to_product
uses to prepare new obs exactly like the hindcast obs, and the bias-correction settings
carried by the hindcast (bias_correction_*, written by bias_correction.bias_correct).
Every variable has 'long_name', 'units' and 'description' attributes (from VARIABLES, and
skill.DESCRIPTIONS for skill_*): values, anomalies and thresholds
are in the units of the variable, percent anomalies in %, probabilities as fractions (0-1),
skill scores as given by skill.hindcast_skill (rmse, bias, crps in the variable's units,
the others dimensionless, units '-').

Later, with the obs of the forecast period aligned with the forecast::

    verify_probabilistic(xr.open_dataset("prob.nc"), obs)     # verify_probabilistic.py
    verify_deterministic(xr.open_dataset("det.nc"), obs)      # verify_deterministic.py
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from toolkit_extension.postprocess import _log
from toolkit_extension.categories import Categories, quiet_nan_warnings
from toolkit_extension.skill import hindcast_skill

# hindcast scores stored by default; roc / bss come per category
PROBABILISTIC_SCORES = ["rpss", "2afc", "roc", "bss", "heidke"]
DETERMINISTIC_SCORES = ["corr", "rmse", "bias"]

HCST_DIM = "hcst_valid_time"


# --------------------------------------------------------------------------- helpers
def _hindcast(da):
    """Hindcast variable with its own time dimension (and no forecast-only coordinates)."""
    return da.drop_vars("target_time", errors="ignore").rename(valid_time=HCST_DIM)


def _skill(scores, **inputs):
    if not scores:
        return {}
    skill = hindcast_skill(list(scores), **inputs)
    return {f"skill_{name}": skill[name] for name in skill}


def _align_inputs(obs, hcst, fcst):
    """
    Forecast lead times the hindcast lacks (the extra lead that staggered initialisations
    can give) are dropped with a message. Otherwise obs, hcst and fcst must have the same
    lead times and grid: an error is raised, since xarray arithmetic would silently drop
    the leads or grid points that differ. Returns fcst.
    """
    f, h = fcst["effective_lead_time"].values, hcst["effective_lead_time"].values
    extra = f[~np.isin(f, h)]
    if extra.size:
        _log(f"products: dropping forecast lead times not in the hindcast: "
             f"{extra.astype('timedelta64[D]').astype(int).tolist()} (days)")
        fcst = fcst.sel(effective_lead_time=f[np.isin(f, h)])
    try:
        xr.align(obs, hcst, fcst, join="exact", exclude=["valid_time", "member"])
    except ValueError as err:
        raise ValueError("obs, hcst and fcst must have the same effective_lead_time, latitude "
                         f"and longitude; align them in the pipeline first\n{err}") from None
    return fcst


def _finish(product, obs, hcst, name, settings, attrs, path):
    """Attributes, tidy coordinates, write."""
    product = _annotate(product, obs.attrs.get("units", ""))
    # regridding and aggregation settings carried by the obs (from align_grid / aggregate),
    # needed to prepare new obs the same way (postprocess.align_obs_to_product)
    product.attrs.update({k: v for k, v in obs.attrs.items()
                          if k.startswith(("regrid_", "aggregation_"))})
    # bias-correction settings carried by the hindcast (from bias_correction.bias_correct),
    # for traceability
    product.attrs.update({k: v for k, v in hcst.attrs.items()
                          if k.startswith("bias_correction_")})
    years = obs["valid_time"].dt.year.values
    product.attrs.update(product=name, units=obs.attrs.get("units", ""),
                         hindcast_period=f"{years.min()}-{years.max()}",
                         n_hindcast_years=int(obs.sizes["valid_time"]), **settings, **(attrs or {}))
    keep = ("category_label", "target_time")
    product = product.drop_vars([c for c in product.coords
                                 if c not in product.dims and c not in keep])
    if path is not None:
        product.to_netcdf(path)
    return product


DATA = "data"           # placeholder below: the units of the variable

# variable: (long_name, units, description); long_name and units label plots.
# skill_* variables get theirs from skill.DESCRIPTIONS.
VARIABLES = {
    "fcst_prob": ("Forecast probability", "-", "forecast probability (0-1) of each category"),
    "fcst_n_members": ("Forecast members", "-", "number of forecast members with data"),
    "hcst_prob": ("Hindcast probability", "-", "hindcast probability (0-1) of each category"),
    "hcst_n_members": ("Hindcast members", "-", "number of hindcast members with data"),
    "obs_cat": ("Observed category", "-",
                "observed category: 1 in the observed category, 0 in the others"),
    "threshold_obs": ("Category threshold (obs)", DATA, "category thresholds for observations"),
    "threshold_model": ("Category threshold (model)", DATA, "category thresholds for the model"),
    "clim_prob": ("Climatological probability", "-",
                  "climatological probability (0-1) of each category"),
    "fcst_mean": ("Forecast ensemble mean", DATA, "forecast ensemble mean"),
    "fcst_anomaly": ("Forecast anomaly", DATA, "forecast ensemble mean - climatology"),
    "fcst_anomaly_pct": ("Forecast anomaly", "%", "100 * forecast anomaly / climatology"),
    "hcst_mean": ("Hindcast ensemble mean", DATA, "hindcast ensemble mean"),
    "hcst_anomaly": ("Hindcast anomaly", DATA, "hindcast ensemble mean - climatology"),
    "hcst_anomaly_pct": ("Hindcast anomaly", "%", "100 * hindcast anomaly / climatology"),
    "obs": ("Observed", DATA, "observations of the hindcast period"),
    "clim_obs": ("Observed climatology", DATA, "mean of the observations over the hindcast years"),
    "clim_model": ("Model climatology", DATA, "mean of the hindcast ensemble mean over the "
                                              "hindcast years"),
}


def _annotate(product, units):
    """long_name, units and description on every variable, replacing inherited attributes
    (e.g. the input data's long_name)."""
    for name, (long_name, u, text) in VARIABLES.items():
        if name in product:
            product[name].attrs = dict(long_name=long_name, units=units if u == DATA else u,
                                       description=text)
    return product


def _drop_member_copies(obs):
    return obs.isel(member=0, drop=True) if "member" in obs.dims else obs


# --------------------------------------------------------------------------- products
@quiet_nan_warnings()
def make_probabilistic_product(obs, hcst, fcst, quantiles=None, *, thresholds=None,
                               edges_from="obs",
                               scores=PROBABILISTIC_SCORES, path=None, attrs=None):
    """
    Category probabilities of forecast and hindcast, and the hindcast skill of the hindcast
    probabilities.

    Parameters
    ----------
    obs, hcst, fcst : xarray.DataArray
        Prepared inputs (see module docstring). They must have the same lead times and
        grid, except that forecast lead times the hindcast lacks are dropped with a
        message.
    quantiles : list of float, optional
        Category boundaries as climatological quantiles, e.g. [1/3, 2/3] (default) or
        [0.1, 0.9].
    thresholds : list of float, optional
        Category boundaries as absolute values instead, e.g. [1.0] (in data units).
    edges_from : {"obs", "model"}, default "obs"
        "obs" (bias-corrected data) or "model" (raw data); see categories.Categories.
    scores : list of str, default PROBABILISTIC_SCORES
        Hindcast scores to store (names from skill.SKILL_SCORES); [] for none.
        Ensemble scores (crpss, groc) can be added too.
    path : str, optional
        If given, the product is also written to this netCDF file.
    attrs : dict, optional
        Extra attributes to record, e.g. {"model": "ECMWF", "variable": "tp"}
        (regridding, aggregation and bias-correction settings are recorded
        automatically).

    Returns
    -------
    xarray.Dataset
        The product (see module docstring).

    Raises
    ------
    ValueError
        If obs, hcst and fcst do not have the same lead times and grid.
    """
    obs = _drop_member_copies(obs)
    fcst = _align_inputs(obs, hcst, fcst)
    cats = Categories(quantiles, edges_from, thresholds=thresholds).fit(obs, hcst)

    hcst_prob, obs_cat = cats.apply(hcst), cats.categorise_obs(obs)
    hcst_n = hcst.count("member")
    skill = _skill(scores, hcst_prob=hcst_prob, obs_cat=obs_cat,
                   clim_prob=cats.clim_prob, n_members=hcst_n, obs=obs, hcst=hcst)

    product = cats.to_dataset()
    product = product.assign(fcst_prob=cats.apply(fcst),
                             fcst_n_members=fcst.count("member"),
                             hcst_prob=_hindcast(hcst_prob),
                             hcst_n_members=_hindcast(hcst_n),
                             obs_cat=_hindcast(obs_cat),
                             obs=_hindcast(obs),
                             **skill)
    return _finish(product, obs, hcst, "probabilistic", {}, attrs, path)   # settings: cats attrs


@quiet_nan_warnings()
def make_deterministic_product(obs, hcst, fcst, *, clim_from="obs",
                               scores=DETERMINISTIC_SCORES, path=None, attrs=None):
    """
    Ensemble means, anomalies and percent anomalies of forecast and hindcast, and the
    hindcast skill of the ensemble mean.

    Parameters
    ----------
    obs, hcst, fcst : xarray.DataArray
        Prepared inputs (see module docstring). They must have the same lead times and
        grid, except that forecast lead times the hindcast lacks are dropped with a
        message.
    clim_from : {"obs", "model"}, default "obs"
        Climatology the anomalies are relative to:

        - "obs": mean of the observations over the hindcast years (bias-corrected data)
        - "model": mean of the hindcast ensemble mean (raw data)

        Future observations are always compared as obs - clim_obs.
    scores : list of str, default DETERMINISTIC_SCORES
        Hindcast scores to store (deterministic family of skill.SKILL_SCORES).
    path, attrs : optional
        As for make_probabilistic_product.

    Returns
    -------
    xarray.Dataset
        The product (see module docstring).

    Raises
    ------
    ValueError
        If clim_from is not "obs" or "model", or if obs, hcst and fcst do not have the
        same lead times and grid.
    """
    if clim_from not in ("obs", "model"):
        raise ValueError("clim_from must be 'obs' or 'model'")
    obs = _drop_member_copies(obs)
    fcst = _align_inputs(obs, hcst, fcst)
    fcst_mean, hcst_mean = fcst.mean("member"), hcst.mean("member")
    clim_obs = obs.mean("valid_time")
    clim = clim_obs if clim_from == "obs" else hcst_mean.mean("valid_time")

    skill = _skill(scores, obs=obs, hcst=hcst)

    def pct(anomaly):
        return 100 * anomaly / clim.where(clim > 0)

    product = xr.Dataset(dict(fcst_mean=fcst_mean, fcst_anomaly=fcst_mean - clim,
                              fcst_anomaly_pct=pct(fcst_mean - clim),
                              hcst_mean=_hindcast(hcst_mean),
                              hcst_anomaly=_hindcast(hcst_mean - clim),
                              hcst_anomaly_pct=_hindcast(pct(hcst_mean - clim)),
                              obs=_hindcast(obs), clim_obs=clim_obs, **skill))
    if clim_from == "model":
        product["clim_model"] = clim
    return _finish(product, obs, hcst, "deterministic", dict(clim_from=clim_from), attrs, path)
