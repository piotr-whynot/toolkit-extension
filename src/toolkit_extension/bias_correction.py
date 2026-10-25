# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

"""
Bias correction of hindcasts and forecasts in effective-lead-time format:

    obs  : (valid_time, effective_lead_time, latitude, longitude)
    hcst : (valid_time, effective_lead_time, member, latitude, longitude)
    fcst : (valid_time, effective_lead_time, member, latitude, longitude)

Internally the data are converted to the corrector format of xsdba_corrector (valid_time
renamed to 'time', masked cells removed) and converted back afterwards.
"""

import inspect
import json
import warnings
from functools import partial

import numpy as np
import xarray as xr

from toolkit_extension.xsdba_corrector import BiasCorrector, XsdbaCorrector, _make_lead_blocks


#this determines with functions will be publicly visible in the installed package
__all__=["bias_correct"]


#######################################################################################
#
# helper functions
#
#######################################################################################

VERBOSE = True

def set_verbose(v):
    """Switch log messages (_log) on or off."""
    global VERBOSE
    VERBOSE = v
    print("logging is {}".format(v))

def _log(msg, force=False):
    """Print `msg` if verbosity is on, or always if force=True."""
    if VERBOSE | force:
        print(msg)



# --------------------------------------------------------------------------- bespoke method
class MeanVarCorrector(BiasCorrector):
    """
    Correction of mean and spread, per lead time and grid point.

    kind="+" (default; temperature and other unbounded variables) - mean/variance scaling:
        corrected = (sim - mean_hist) / std_hist * std_ref + mean_ref
    kind="*" (precipitation and other non-negative variables) - power transformation
    (Leander & Buishand 2007, J. Hydrol. 332):
        corrected = a * sim**b
        b is chosen so that the coefficient of variation (std/mean) of hist**b equals that of
        ref, then a so that the means match. Output is never negative and zeros stay zero.
        b is searched in [b_min, b_max]; where the target is outside that range, b is clipped
        and the spread is only partly corrected. Negative input values are set to 0.

    Statistics always pool all years (and, for the hindcast, all members). Pooling across
    lead times, as in XsdbaCorrector:
        "lead"    : statistics per lead time
        "rolling" : per lead time, from a centred window of `window_size` leads
                    (truncated at the first and last leads)
        "blocks"  : one set of statistics per block of `window_size` leads, used for every
                    lead in the block; a short last block is merged into the previous one

    Parameters
    ----------
    pooling : "lead" (default), "rolling" or "blocks", as above
    window_size : leads per window or block (default 7); ignored for pooling="lead"
    kind : "+" (default) or "*", as above
    b_min, b_max : search range of the exponent b for kind="*" (default 0.1 to 4.0)
    """

    POOLINGS = ("lead", "rolling", "blocks")

    def __init__(self, pooling="lead", window_size=7, kind="+", b_min=0.1, b_max=4.0):
        if pooling not in self.POOLINGS:
            raise ValueError(f"pooling must be one of {self.POOLINGS}")
        if kind not in ("+", "*"):
            raise ValueError("kind must be '+' or '*'")
        self.pooling, self.window_size, self.kind = pooling, int(window_size), kind
        self.b_grid = np.linspace(b_min, b_max, int(round((b_max - b_min) / 0.05)) + 1)

    def _stats(self, da, dims, power=None):
        """Mean and std of da (or da**power) over `dims`, pooled across leads per `pooling`.
        `power` may be a scalar or an array along (lead, ...); in a rolling window each
        window uses the power of its central lead."""
        L = "effective_lead_time"
        if self.pooling == "lead":
            x = da if power is None else da ** power
            return x.mean(dims), x.std(dims)
        if self.pooling == "rolling":
            w = da.rolling({L: self.window_size}, center=True).construct("_win")
            w = w if power is None else w ** power
            return w.mean(dims + ["_win"]), w.std(dims + ["_win"])
        # blocks: statistics per block, then mapped back onto every lead of the block
        # (power is constant within a block, because the statistics it is fitted to are)
        x = da if power is None else da ** power
        blocks = _make_lead_blocks(da.sizes[L], self.window_size)
        block_of_lead = np.concatenate([np.full(len(b), i) for i, b in enumerate(blocks)])
        idx = xr.DataArray(block_of_lead, dims=L)
        g = x.assign_coords(_block=(L, block_of_lead)).groupby("_block")
        mean = g.mean(dims + [L]).isel(_block=idx).drop_vars("_block")
        std = g.std(dims + [L]).isel(_block=idx).drop_vars("_block")
        return mean.assign_coords({L: da[L]}), std.assign_coords({L: da[L]})

    def _fit_power(self, ref, hist, hdims):
        """Exponent b with CV(hist**b) = CV(ref), by linear interpolation on b_grid."""
        r_mean, r_std = self._stats(ref, ["time"])
        target = r_std / r_mean
        cv = xr.concat([(lambda m, s: s / m)(*self._stats(hist, hdims, power=b))
                        for b in self.b_grid], dim="_b")          # increases with b
        target = target.broadcast_like(cv.isel(_b=0))
        cvv = cv.transpose("_b", *target.dims).values
        tv = target.values
        above = cvv >= tv[None]
        k = np.clip(np.argmax(above, axis=0), 1, len(self.b_grid) - 1)  # first b with cv >= target
        k = np.where(above.any(axis=0), k, len(self.b_grid) - 1)
        lo = np.take_along_axis(cvv, (k - 1)[None], 0)[0]
        hi = np.take_along_axis(cvv, k[None], 0)[0]
        frac = np.clip(np.where(hi > lo, (tv - lo) / (hi - lo), 1.0), 0, 1)
        b = self.b_grid[k - 1] + frac * (self.b_grid[k] - self.b_grid[k - 1])
        b = np.where(np.isfinite(tv) & np.isfinite(cvv).all(axis=0), b, np.nan)
        b = xr.DataArray(b, dims=target.dims, coords=target.coords)
        h_mean_b, _ = self._stats(hist, hdims, power=b)
        a = (r_mean / h_mean_b).where(h_mean_b > 0, 0.0)
        return a, b

    def fit(self, ref, hist):
        self._lead = hist["effective_lead_time"].values
        hdims = ["time", "member"]
        if self.kind == "*":
            hist, ref = hist.clip(min=0), ref.clip(min=0)
            self.a, self.b = self._fit_power(ref, hist, hdims)
        else:
            self.h_mean, self.h_std = self._stats(hist, hdims)
            self.r_mean, self.r_std = self._stats(ref, ["time"])
        return self

    def apply(self, sim):
        # arithmetic aligns on lead labels; mismatched leads would silently disappear
        if not np.array_equal(sim["effective_lead_time"].values, self._lead):
            raise ValueError("'effective_lead_time' of sim differs from the training data")
        if self.kind == "*":
            out = self.a * sim.clip(min=0) ** self.b
        else:
            out = (sim - self.h_mean) / self.h_std * self.r_std + self.r_mean
        return out.transpose(*sim.dims).assign_attrs(sim.attrs)




# --------------------------------------------------------------------------- registry
REGISTRY: dict[str, object] = {}


def register(name, corrector):
    """Make a corrector available to bias_correct under `name`.
    corrector : class (or callable) returning a BiasCorrector; called with the user's
                options (pooling, window_size and method options)."""
    REGISTRY[name] = corrector


for _m in ("eqm", "scaling", "dqm"):
    register(_m, partial(XsdbaCorrector, _m))
register("meanvar", MeanVarCorrector)



def _match_leads(fcst, hcst):
    """Restrict the forecast to the hindcast's lead times. Staggered inits can leave the
    forecast with extra leads; leads the hindcast has but the forecast lacks are an error."""
    f, h = fcst["effective_lead_time"].values, hcst["effective_lead_time"].values
    missing = h[~np.isin(h, f)]
    if missing.size:
        raise ValueError(f"forecast lacks hindcast lead times {missing}")
    extra = f[~np.isin(f, h)]
    if extra.size:
        _log(f"bias_correct: dropping {extra.size} forecast lead times not in the hindcast: {extra}")
        fcst = fcst.sel(effective_lead_time=h)
    return fcst


def _to_points(da, mask):
    """
    Replace (latitude, longitude) with a single 'point' dimension holding only cells where
    `mask` is True (e.g. land). latitude/longitude are kept as plain coordinates along
    'point' (no MultiIndex, which xsdba cannot handle). Skipping masked cells avoids
    all-NaN warnings and wasted work.
    """
    st = da.stack(point=("latitude", "longitude"))
    keep = mask.stack(point=("latitude", "longitude")).values.astype(bool)
    lats, lons = st["latitude"].values[keep], st["longitude"].values[keep]
    out = st.drop_vars(["point", "latitude", "longitude"]).isel(point=np.flatnonzero(keep))
    return out.assign_coords(latitude=("point", lats), longitude=("point", lons))


def _from_points(da, grid):
    """
    Inverse of _to_points: put 'point' data back onto the full latitude/longitude grid,
    NaN in masked cells.

    grid : any array with latitude and longitude dimensions that defines the full grid,
           e.g. the mask. The other dimensions keep their order; (latitude, longitude)
           take the place of 'point'.
    """
    order = []
    for d in da.dims:
        order.extend(["latitude", "longitude"] if d == "point" else [d])
    out = da.set_index(point=["latitude", "longitude"]).unstack("point")
    out = out.reindex(latitude=grid["latitude"].values, longitude=grid["longitude"].values)
    return out.transpose(*order).assign_coords(latitude=grid["latitude"],
                                               longitude=grid["longitude"])



# --------------------------------------------------------------------------- format conversion
def _to_bc(da, mask, units):
    """Package format -> corrector format."""
    # non-index coords and member labels go (member labels would make xsdba align the
    # forecast against the hindcast's members); they are restored from the original later
    out = da.reset_coords(drop=True).drop_vars("member", errors="ignore")
    out = out.rename(valid_time="time")
    if units is not None:
        out = out.assign_attrs(units=units)
    if mask is None:
        return out
    return _to_points(out, mask)


def _from_bc(da, like, mask):
    """Corrector format -> package format of `like` (the original array)."""
    if mask is not None:
        da = _from_points(da, grid=mask)
    out = da.rename(time="valid_time")
    out = out.assign_coords({c: like[c] for c in like.coords if c not in out.coords})
    out = out.transpose(*like.dims)
    if out.sizes != like.sizes:
        raise RuntimeError(f"size mismatch: {dict(out.sizes)} vs {dict(like.sizes)}")
    return out


def _check_units(obs, hcst, units):
    """With units=None both inputs need a 'units' attribute, and the same one."""
    if units is not None:
        return
    u_obs, u_hc = obs.attrs.get("units"), hcst.attrs.get("units")
    if u_obs is None or u_hc is None:
        raise ValueError("obs and hcst need a 'units' attribute; pass units=...")
    if u_obs != u_hc:
        raise ValueError(f"obs units '{u_obs}' differ from hcst units '{u_hc}'; convert the "
                         "data or set the units when reading them (read_netcdf(..., units=...))")


def _accepted_options(factory, options):
    """Split options into those the corrector accepts and those it doesn't."""
    params = inspect.signature(factory).parameters
    if any(p.kind is p.VAR_KEYWORD for p in params.values()):    # takes **kwargs: pass all
        return options, {}
    used = {k: v for k, v in options.items() if k in params}
    unused = {k: v for k, v in options.items() if k not in params}
    return used, unused



# --------------------------------------------------------------------------- entry point
def bias_correct(obs, hcst, fcst=None, method="eqm", *, cv=True, pooling="lead", window_size=7,
                 leave_out=1, mask="auto", units=None, return_model=False, **options):
    """
    Bias-correct a hindcast and a forecast.

    Parameters
    ----------
    obs   : (valid_time, effective_lead_time, latitude, longitude), aligned to hcst
            (align_obs_to_model_time); a member dimension of identical copies is allowed
    hcst  : (valid_time, effective_lead_time, member, latitude, longitude)
    fcst  : (valid_time, effective_lead_time, member, latitude, longitude), optional; any
            number of members; lead times the hindcast lacks are dropped
    method : name in REGISTRY ("eqm", "scaling", "dqm", "meanvar", or one added with register)
    cv : if True (default) the returned hindcast is cross-validated (`leave_out` years at a
         time), as needed for skill assessment; if False it is corrected in-sample. The
         forecast is always corrected with the model fitted on all years.
    pooling : "lead" (default), "rolling" or "blocks": how data are pooled across lead times
         (see XsdbaCorrector / MeanVarCorrector)
    window_size : leads per window or block for pooling="rolling" / "blocks" (default 7)
    leave_out : number of consecutive years left out at a time when cv=True (default 1)
    mask : "auto" (cells where obs have any data), a boolean (latitude, longitude) array, or
           None. Cells outside the mask are not corrected and are NaN in the output.
    units : replaces the 'units' attribute of obs, hcst and fcst (required by the xsdba
           methods). If None, obs and hcst must have the same 'units' attribute.
    return_model : also return the corrector fitted on all years.
    **options : method options, e.g. nquantiles, kind, interp (see XsdbaCorrector and
           MeanVarCorrector). Options the method does not accept are ignored with a warning.

    Returns
    -------
    hcst_bc, fcst_bc            (fcst_bc is None if fcst is None)
    hcst_bc, fcst_bc, model     if return_model
    The settings are recorded as attributes of hcst_bc and fcst_bc: bias_correction_method,
    bias_correction_pooling, bias_correction_window, bias_correction_cv,
    bias_correction_leave_out (hindcast, if cv) and bias_correction_options (JSON).
    """
    if method not in REGISTRY:
        raise ValueError(f"unknown method '{method}'; available: {list(REGISTRY)}")

    if "member" in obs.dims:                                    # broadcast copies -> drop
        obs = obs.isel(member=0, drop=True)
    if isinstance(mask, str) and mask == "auto":
        mask = obs.notnull().any(["valid_time", "effective_lead_time"])
    _check_units(obs, hcst, units)

    ref = _to_bc(obs, mask, units)
    hist = _to_bc(hcst, mask, units)

    factory = REGISTRY[method]
    options = dict(options, pooling=pooling, window_size=window_size)
    used, unused = _accepted_options(factory, options)
    if unused:
        warnings.warn(f"method '{method}' ignores options {sorted(unused)}")
    model = factory(**used)

    # the model fitted on all years: for the in-sample hindcast, the forecast, or to return
    if not cv or fcst is not None or return_model:
        model.fit(ref, hist)

    if cv:
        hist_bc = model.fit_apply_loo(ref, hist, leave_out)
    else:
        hist_bc = model.apply(hist)
    hcst_bc = _from_bc(hist_bc, hcst, mask)

    # settings, for traceability
    method_options = {k: v for k, v in used.items() if k not in ("pooling", "window_size")}
    settings = dict(bias_correction_method=method,
                    bias_correction_pooling=pooling,
                    bias_correction_window=int(window_size),
                    bias_correction_options=json.dumps(method_options, default=str))
    hcst_bc = hcst_bc.assign_attrs(settings, bias_correction_cv=str(cv))
    if cv:
        hcst_bc = hcst_bc.assign_attrs(bias_correction_leave_out=int(leave_out))

    fcst_bc = None
    if fcst is not None:
        fcst = _match_leads(fcst, hcst)
        sim = _to_bc(fcst, mask, units)
        fcst_bc = _from_bc(model.apply(sim), fcst, mask)
        fcst_bc = fcst_bc.assign_attrs(settings, bias_correction_cv="False")

    if return_model:
        return hcst_bc, fcst_bc, model
    return hcst_bc, fcst_bc
