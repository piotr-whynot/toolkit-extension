"""
Verification of an issued probabilistic (category) forecast against the observations that
have since become available, following the WMO guidance on verifying operational seasonal
forecasts (WMO-No. 1220).

    product = xr.open_dataset("prob.nc")             # from make_probabilistic_product
    skill = verify_probabilistic(product, obs)
    skill = verify_probabilistic(product, obs, ["heidke_hits", "heidke"])

Inputs, alignment rules, area weighting: see verification_common.py.
The obs are categorised with the thresholds stored in the product, exactly as in the
hindcast; the reference forecast is the stored clim_prob (quantile spacing, or the observed
frequency for absolute thresholds).

Units: obs in the units of the variable; obs_percentile, interest_rate and
effective_interest_rate in %; ignorance in bits; probabilities as fractions (0-1); all other
scores dimensionless (units '-'). Every output variable carries 'units'.

Maps (valid_time, lead, latitude, longitude)
    obs             the observations of the forecast period, as verified against
    obs_percentile  the observation as a percentile (0-100) of the hindcast-period obs
                    stored in the product, at that grid point and lead time (hindcast values
                    equal to the obs count half)
    obs_category    index of the observed category
    prob_observed   probability the forecast gave to the observed category
    heidke_hits     1 if the most probable category was observed, 0 if not; with ties
                    1/n for each of the n tied categories (e.g. 0.5 for two)
    ignorance       -log2(probability of the observed category), in bits; lower is better
                    (climatology for terciles: log2(3) = 1.58)
    interest_rate   100 * (p_obs / clim_obs - 1), in %: the return on a bet placed on the
                    forecast probabilities against odds set by climatology (0 = no gain,
                    100 = money doubled, -100 = all lost)
    rps             fair ranked probability score (lower is better)

Scores over the map (valid_time, lead[, category])
    heidke          Heidke skill score of the most probable category: (H - E) / (1 - E),
                    H = area fraction of hits (heidke_hits), E = hits expected by chance =
                    climatological probability of the forecast category (1/3 for terciles);
                    0 = chance, 1 = perfect, negative = worse than chance
    mean_ignorance  area mean of ignorance
    effective_interest_rate
                    the average return over the map (geometric mean, as compounded
                    interest), in % = 100 * (2^(mean ignorance of climatology - mean
                    ignorance) - 1)
    rpss            1 - mean(RPS) / mean(RPS of climatology)
    bss             per category: 1 - mean(BS) / mean(BS of climatology) of the event
                    "obs in this category" (fair Brier score)
With two categories, rpss equals bss of either category.

Ignorance and interest rate are infinite / -100% where the forecast gave 0 probability to
what happened, which ensemble fractions often do (e.g. no member above the threshold). For
these two scores the forecast probabilities are therefore those of the ensemble plus one
member drawn from climatology: (M * p + clim_prob) / (M + 1), M = number of members.
Where the observed category was never observed in the hindcast years (clim_prob 0, possible
with absolute thresholds) there is no climatological reference: ignorance and interest rate
are NaN there, and those grid points are left out of mean_ignorance and
effective_interest_rate. rpss and bss are NaN where their reference score is 0, heidke where
the chance of a hit is 1 (climatology certain of one category everywhere).
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from toolkit_extension.categories import Categories, quiet_nan_warnings
from toolkit_extension.verification_common import (DATA_UNITS, area_mean, check_product, finish, obs_percentile, plain,
                                 hindcast_obs, prepare_obs, select_scores)

# name: (long_name, units, description); long_name and units label plots
PROBABILISTIC_SCORES = {
    "obs": ("Observed", DATA_UNITS, "observations of the forecast period (the obs verified "
                                    "against)"),
    "obs_percentile": ("Observed percentile", "%",
                       "observation as a percentile (0-100) of the hindcast-period observations"),
    "obs_category": ("Observed category", "-", "index of the observed category"),
    "prob_observed": ("Probability of observed category", "-",
                      "forecast probability (0-1) of the observed category"),
    "heidke_hits": ("Hits of the most probable category", "-",
                    "1 = most probable category observed, 0 = not, 1/n for n tied categories"),
    "heidke": ("Heidke skill score", "-", "Heidke skill score of the most probable category "
                                          "over the map, chance = climatological probability"),
    "ignorance": ("Ignorance", "bit", "-log2(probability of the observed category), forecast with "
                                      "one climatological member added (lower is better)"),
    "mean_ignorance": ("Mean ignorance", "bit", "area-weighted mean ignorance"),
    "interest_rate": ("Interest rate", "%", "100 * (probability of the observed category, forecast "
                                            "with one climatological member added / its "
                                            "climatological probability - 1)"),
    "effective_interest_rate": ("Effective interest rate", "%",
                                "100 * (2^(mean ignorance of climatology - mean ignorance) - 1)"),
    "rps": ("Ranked probability score", "-", "fair ranked probability score (lower is better)"),
    "rpss": ("Ranked probability skill score", "-",
             "1 - mean(RPS) / mean(RPS of climatology), area-weighted"),
    "bss": ("Brier skill score", "-",
            "1 - mean(BS) / mean(BS of climatology) per category, area-weighted"),
}


@quiet_nan_warnings()
def verify_probabilistic(product, obs, scores=None):
    """
    Skill of the probabilistic forecast in `product` (see module docstring).

    product : Dataset from make_probabilistic_product (or the saved file)
    obs     : (valid_time, lead, latitude, longitude), aligned with the forecast, e.g. by
              postprocess.align_obs_to_product; extra trailing lead times are dropped
    scores  : names from PROBABILISTIC_SCORES; None = all

    Returns a Dataset with one variable per score.
    """
    check_product(product, "probabilistic")
    names = select_scores(scores, PROBABILISTIC_SCORES)
    obs = prepare_obs(obs, product, "fcst_prob")
    hcst_obs = hindcast_obs(product)                 # hindcast-period obs, for percentiles

    p = plain(product["fcst_prob"])
    o = plain(Categories.from_dataset(product).categorise_obs(obs))
    c = (xr.zeros_like(p) + plain(product["clim_prob"])).where(p.notnull())
    m = plain(product["fcst_n_members"])
    valid = p.notnull().all("category") & o.notnull().all("category")

    def rps(prob, n=None):
        """RPS per grid point; fair (ensemble-size corrected) if n is given."""
        F, O = prob.cumsum("category"), o.cumsum("category")
        r = ((F - O) ** 2).sum("category")
        if n is not None:
            r = r - (F * (1 - F)).sum("category") / (n - 1)
        return r.where(valid)

    # most probable category, ties shared
    top = abs(p - p.max("category")) < 1e-9
    n_top = top.sum("category")
    hits = ((top * o).sum("category") / n_top).where(valid)
    chance = ((top * c).sum("category") / n_top).where(valid)

    # probability of the observed category: forecast (plus one climatological member) and
    # climatology, for ignorance and interest rate. Where the observed category was never
    # observed in the hindcast years (clim_prob 0, possible with absolute thresholds), there
    # is no climatological reference and the scores would be infinite: NaN there.
    c_obs = (c * o).sum("category").where(valid)
    seen = valid & (c_obs > 0)
    p_obs = (((m * p + c) / (m + 1)) * o).sum("category").where(seen)
    c_obs = c_obs.where(seen)

    out = {}
    for name in names:
        if name == "obs":
            out[name] = obs
        elif name == "obs_percentile":
            out[name] = obs_percentile(obs, hcst_obs)
        elif name == "obs_category":
            out[name] = o.fillna(-1).argmax("category").where(valid)
        elif name == "prob_observed":
            out[name] = (p * o).sum("category").where(valid)
        elif name == "heidke_hits":
            out[name] = hits
        elif name == "heidke":
            H, E = area_mean(hits, valid), area_mean(chance, valid)
            out[name] = ((H - E) / (1 - E)).where(E < 1)   # E = 1: chance always hits
        elif name == "ignorance":
            out[name] = -np.log2(p_obs)
        elif name == "mean_ignorance":
            out[name] = area_mean(-np.log2(p_obs), seen)
        elif name == "interest_rate":
            out[name] = 100 * (p_obs / c_obs - 1)
        elif name == "effective_interest_rate":
            out[name] = 100 * (2 ** area_mean(np.log2(p_obs) - np.log2(c_obs), seen) - 1)
        elif name == "rps":
            out[name] = rps(p, m)
        elif name == "rpss":
            ref = area_mean(rps(c), valid)
            out[name] = (1 - area_mean(rps(p, m), valid) / ref).where(ref > 0)
        elif name == "bss":
            bs = (p - o) ** 2 - p * (1 - p) / (m - 1)
            ref = area_mean((c - o) ** 2, valid)          # 0: climatology certain and right
            out[name] = (1 - area_mean(bs, valid) / ref).where(ref > 0)
    return finish(out, names, PROBABILISTIC_SCORES, product)
