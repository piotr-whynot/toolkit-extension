"""
Shared parts of forecast verification (verify_deterministic.py, verify_probabilistic.py):
input checks, area-weighted means over the map, the observed percentile.

Inputs of both verification functions
    product : the saved product (make_deterministic_product / make_probabilistic_product).
              Its hindcast-period obs ('obs', along hcst_valid_time) are the climatology for
              the observed percentile.
    obs     : (valid_time, effective_lead_time, latitude, longitude), aligned with the
              forecast in the pipeline (same valid_time, lead times, grid, aggregation and
              units). Trailing lead times the product lacks (the extra lead that the
              forecast's initialisation lag can give) are dropped; anything else that does
              not match raises an error.

A single forecast has one value per grid point, so results are either
    maps    (valid_time, effective_lead_time, latitude, longitude)
    scores  (valid_time, effective_lead_time): over the map, weighted by cos(latitude),
            using only grid points with both forecast and obs
"""

from __future__ import annotations

import numpy as np
import xarray as xr

SPACE = ["latitude", "longitude"]


def _log(msg):          # replace with your package's _log
    print(msg)


# --------------------------------------------------------------------------- input checks
def plain(da):
    """Without non-index coordinates (e.g. target_time), except category_label."""
    return da.drop_vars([c for c in da.coords if c not in da.dims and c != "category_label"])


def check_product(product, kind):
    """Raise unless `product` is a product of the given kind."""
    found = product.attrs.get("product")
    if found != kind:
        raise ValueError(f"expected a {kind} product, got '{found}' "
                         f"(use verify_{found} for it)" if found else
                         f"expected a {kind} product (attribute 'product' missing)")


def select_scores(scores, available):
    """The requested score names (all if None), checked against `available`."""
    names = list(available) if scores is None else list(scores)
    unknown = [s for s in names if s not in available]
    if unknown:
        raise ValueError(f"unknown scores {unknown}; available: {list(available)}")
    return names


def _check_aligned(da, ref, what):
    try:
        xr.align(da, ref, join="exact", exclude=["category"])
    except ValueError as err:
        raise ValueError(f"{what}; align them in the pipeline before verification\n{err}") from None


def prepare_obs(obs, product, fcst_var):
    """
    Check that obs are aligned with the forecast `product[fcst_var]` (valid_time, lead
    times, grid, units). The only adjustment: trailing lead times the product lacks are
    dropped.
    """
    if "member" in obs.dims:                                      # broadcast copies
        obs = obs.isel(member=0, drop=True)
    leads = product["effective_lead_time"].values
    n = len(leads)
    if (obs.sizes["effective_lead_time"] > n
            and np.array_equal(obs["effective_lead_time"].values[:n], leads)):
        extra = obs["effective_lead_time"].values[n:]
        _log(f"verification: dropping obs lead times not in the product: "
             f"{extra.astype('timedelta64[D]').astype(int).tolist()} (days)")
        obs = obs.isel(effective_lead_time=slice(0, n))
    _check_aligned(obs, product[fcst_var],
                   "obs are not aligned with the forecast (valid_time, lead times or grid)")
    u_obs, u_fc = obs.attrs.get("units"), product.attrs.get("units")
    if u_obs and u_fc and u_obs != u_fc:
        raise ValueError(f"obs units '{u_obs}' differ from the product's '{u_fc}'")
    return plain(obs)


def hindcast_obs(product):
    """The hindcast-period obs stored in the product, (hcst_valid_time, lead, lat, lon)."""
    if "obs" not in product:
        raise ValueError("the product has no hindcast-period obs ('obs'); it was made by an "
                         "older version of make_*_product - regenerate it")
    return plain(product["obs"])


# --------------------------------------------------------------------------- calculations
def area_mean(x, valid):
    """cos(latitude)-weighted mean over the map, over valid points only."""
    w = np.cos(np.deg2rad(x["latitude"]))
    return x.where(valid).weighted(w).mean(SPACE)


def area_weights(da):
    """cos(latitude) weights on the (latitude, longitude) grid of `da`."""
    return np.cos(np.deg2rad(da["latitude"])) * xr.ones_like(da["longitude"], dtype=float)


def obs_percentile(obs, clim):
    """obs as a percentile (0-100) of the hindcast obs `clim` along hcst_valid_time, per
    grid point and lead time; hindcast values equal to the obs count half."""
    below = (clim < obs).sum("hcst_valid_time")
    equal = (clim == obs).sum("hcst_valid_time")
    n = clim.count("hcst_valid_time")
    return (100 * (below + 0.5 * equal) / n).where((n > 0) & obs.notnull())


DATA_UNITS = "data"      # placeholder in the score tables: the units of the variable


def finish(results, names, table, product):
    """Dataset of the requested results, in the requested order, with long_name, units and
    description from `table` {name: (long_name, units, description)}; DATA_UNITS is replaced
    by the product's units. Inherited attributes are dropped (they would mislabel plots)."""
    out = {}
    for n in names:
        long_name, units, text = table[n]
        if units == DATA_UNITS:
            units = product.attrs.get("units", "")
        da = results[n].rename(n)
        da.attrs = dict(long_name=long_name, units=units, description=text)
        out[n] = da
    return xr.Dataset(out).assign_attrs(product=product.attrs.get("product", ""),
                                        units=product.attrs.get("units", ""))
