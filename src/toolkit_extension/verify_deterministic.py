"""
Verification of an issued deterministic forecast (ensemble mean) against the observations
that have since become available::

    product = xr.open_dataset("det.nc")              # from make_deterministic_product
    skill = verify_deterministic(product, obs)
    skill = verify_deterministic(product, obs, ["acc", "rmse"])

Inputs, alignment rules, area weighting: see verification_common.py.
Anomalies are relative to the hindcast climatology stored in the product (clim_obs for the
observations; the product's own fcst_anomaly for the forecast).

Units: obs, error, obs_anomaly, bias, mae and rmse are in the units of the variable;
obs_percentile, error_pct and bias_pct in %; acc, msss and sign_agreement are dimensionless
(units '-'). Every output variable carries 'units'.

Maps (valid_time, lead, latitude, longitude)::

    obs             the observations of the forecast period, as verified against
    obs_percentile  the observation as a percentile (0-100) of the hindcast-period obs
                    stored in the product, at that grid point and lead time (hindcast values
                    equal to the obs count half)
    error           fcst_mean - obs
    error_pct       100 * (fcst_mean - obs) / clim_obs: the error as a percentage of the
                    climatology (= forecast minus observed percent anomaly). Relative to the
                    climatology rather than to the observed value, which can be 0 (e.g. a
                    dry week); NaN where the climatology is 0 or negative
    obs_anomaly     obs - clim_obs, to plot next to fcst_anomaly

Scores over the map (valid_time, lead)::

    bias            mean error
    bias_pct        percent bias: 100 * mean error / mean obs (NaN if mean obs <= 0)
    mae             mean absolute error
    rmse            root-mean-square error
    acc             anomaly correlation: correlation over the map between forecast and
                    observed anomalies (centred, i.e. spatial means removed)
    msss            mean-squared skill score of the anomalies against climatology (zero
                    anomaly): 1 - MSE(fcst anomaly) / mean(obs anomaly^2);
                    0 = no better than climatology, 1 = perfect
    sign_agreement  fraction of the area where forecast and observed anomalies have the same
                    sign (0.5 = chance for a symmetric variable)

Note: forecast and observed anomalies share the error of the climatology estimated from
the hindcast years, which correlates them slightly; with 20 years an unrelated forecast
gets ACC of about +0.07 on average rather than 0.
"""

from __future__ import annotations

import numpy as np
import xskillscore as xs

from toolkit_extension.categories import quiet_nan_warnings
from toolkit_extension.verification_common import (DATA_UNITS, SPACE, area_mean, area_weights, check_product, finish,
                                 obs_percentile, plain, hindcast_obs, prepare_obs,
                                 select_scores)

# name: (long_name, units, description); long_name and units label plots;
# DATA_UNITS = the units of the variable (product attribute)
DETERMINISTIC_SCORES = {
    "obs": ("Observed", DATA_UNITS, "observations of the forecast period (the obs verified "
                                    "against)"),
    "obs_percentile": ("Observed percentile", "%",
                       "observation as a percentile (0-100) of the hindcast-period observations"),
    "error": ("Forecast error", DATA_UNITS, "fcst_mean - obs"),
    "error_pct": ("Forecast error", "%", "100 * (fcst_mean - obs) / clim_obs"),
    "obs_anomaly": ("Observed anomaly", DATA_UNITS, "obs - clim_obs"),
    "bias": ("Bias", DATA_UNITS, "area-weighted mean of fcst_mean - obs"),
    "bias_pct": ("Percent bias", "%", "100 * area-weighted mean of (fcst_mean - obs) / "
                                      "area-weighted mean of obs"),
    "mae": ("Mean absolute error", DATA_UNITS, "area-weighted mean absolute error"),
    "rmse": ("RMSE", DATA_UNITS, "area-weighted root-mean-square error"),
    "acc": ("Anomaly correlation", "-",
            "area-weighted correlation of forecast and observed anomalies over the map"),
    "msss": ("Mean-squared skill score", "-",
             "1 - MSE(fcst anomaly) / MSE(zero anomaly), area-weighted (1 = perfect, "
             "0 = climatology)"),
    "sign_agreement": ("Anomaly sign agreement", "-",
                       "area fraction (0-1) where forecast and observed anomalies have the "
                       "same sign"),
}


@quiet_nan_warnings()
def verify_deterministic(product, obs, scores=None):
    """
    Skill of the deterministic forecast in `product` (see module docstring).

    Parameters
    ----------
    product : xarray.Dataset
        Product from make_deterministic_product (or the saved file).
    obs : xarray.DataArray
        (valid_time, lead, latitude, longitude), aligned with the forecast, e.g. by
        postprocess.align_obs_to_product; extra trailing lead times are dropped.
    scores : list of str, optional
        Names from DETERMINISTIC_SCORES; None (default) = all.

    Returns
    -------
    xarray.Dataset
        One variable per score.
    """
    check_product(product, "deterministic")
    names = select_scores(scores, DETERMINISTIC_SCORES)
    obs = prepare_obs(obs, product, "fcst_mean")
    hcst_obs = hindcast_obs(product)                 # hindcast-period obs, for percentiles
    clim = plain(product["clim_obs"])                # their mean, for anomalies

    fc, fa = plain(product["fcst_mean"]), plain(product["fcst_anomaly"])
    oa = obs - clim
    valid = fc.notnull() & obs.notnull()
    err = fc - obs

    out = {}
    for name in names:
        if name == "obs":
            out[name] = obs
        elif name == "obs_percentile":
            out[name] = obs_percentile(obs, hcst_obs)
        elif name == "error":
            out[name] = err
        elif name == "error_pct":
            out[name] = 100 * err / clim.where(clim > 0)
        elif name == "obs_anomaly":
            out[name] = oa
        elif name == "bias":
            out[name] = area_mean(err, valid)
        elif name == "bias_pct":
            mean_obs = area_mean(obs, valid)
            out[name] = 100 * area_mean(err, valid) / mean_obs.where(mean_obs > 0)
        elif name == "mae":
            out[name] = area_mean(abs(err), valid)
        elif name == "rmse":
            out[name] = np.sqrt(area_mean(err ** 2, valid))
        elif name == "acc":
            out[name] = xs.pearson_r(oa.where(valid), fa.where(valid), dim=SPACE,
                                     weights=area_weights(fc), skipna=True)
        elif name == "msss":
            out[name] = 1 - area_mean((fa - oa) ** 2, valid) / area_mean(oa ** 2, valid)
        elif name == "sign_agreement":
            out[name] = area_mean((np.sign(fa) == np.sign(oa)).astype(float), valid)
    return finish(out, names, DETERMINISTIC_SCORES, product)
