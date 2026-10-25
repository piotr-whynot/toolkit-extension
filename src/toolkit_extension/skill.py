"""
Hindcast skill, organised by score family.

Step 1 - prepare the inputs a family needs (once, reused by all its scores):

    ev  = make_event(obs, hcst, quantile=0.9, side="above")       # e.g. > 90th percentile
    cat = make_categories(obs, hcst, quantiles=[1/3, 2/3])        # e.g. terciles
    ens = dict(obs=obs, hcst=hcst)                                # full ensemble

Step 2 - compute scores; the score name decides the family (SKILL_SCORES):

    hindcast_skill("bss", **ev)
    hindcast_skill(["rps", "rpss", "heidke"], **cat)               # list -> Dataset
    hindcast_skill(["crpss", "groc", "corr"], **ens)

Families and their inputs
-------------------------
prob_binary   : hcst_prob, obs_cat (0/1), clim_prob, n_members   <- make_event
                or the output of make_categories: one score per category
                (e.g. ROC of the below-, near- and above-normal tercile)
categorical   : hcst_prob, obs_cat (one-hot), clim_prob, n_members, with dim "category"
                                                                   <- make_categories
ensemble_obs  : hcst (members), obs
deterministic : hcst (members; the ensemble mean is used), obs
Inputs a family does not need are ignored, so all inputs can be passed at once
(as the functions in products.py do).

All scores are computed over valid_time (the hindcast years) for every lead time and grid
point.
Each result carries 'long_name', 'units' and 'description' attributes from DESCRIPTIONS:
crps, rmse and bias are in the units of the variable (taken from obs), all other scores are
dimensionless (units '-').
Thresholds come from categories.Categories, which also makes the forecast products, so
hindcast skill and products always use the same thresholds.

Conventions
-----------
- Thresholds (event / category edges) are absolute values, or quantiles of the observations
  over all years, per lead time and grid point; the hindcast is categorised with the same thresholds
  (edges_from="obs", for bias-corrected data) or with quantiles of its own climatology,
  pooled over years and members (edges_from="model", for raw data). Thresholds are fixed,
  not leave-one-out: a leave-one-out threshold shared by forecast and observed event depends
  on the observation being verified and gives spurious skill (ROC ~0.62 for a no-skill
  forecast in tests).
- Reference forecast for BSS / RPSS: the climatological probabilities (1 - q for "> q-th
  percentile", the quantile spacing for categories). For CRPSS: the observations of all other
  years as an ensemble (leave-one-year-out).
  The climatological probabilities match the observed frequencies exactly (the thresholds
  come from the same years), whereas the hindcast's are subject to threshold sampling error,
  so a no-skill hindcast gets BSS / RPSS of about -1/n_years (-0.03 for 30 years), not 0.
- Brier, RPS and CRPS are "fair" (Ferro 2014, QJRMS 140): corrected for ensemble size, so that
  ensembles of different sizes, and the climatological reference, are comparable.
- NaN members (e.g. staggering gaps) are ignored when computing probabilities.

Adding a score: add its name to the family list in SKILL_SCORES, a line to DESCRIPTIONS and a
branch in that family's block of hindcast_skill.
"""

from __future__ import annotations

import numpy as np
import xarray as xr
import xskillscore as xs

from toolkit_extension.categories import Categories, order_dims, quiet_nan_warnings


SKILL_SCORES = {
    "prob_binary": ["bs", "bss", "roc"],
    "categorical": ["rps", "rpss", "heidke", "2afc"],
    "ensemble_obs": ["crps", "crpss", "groc"],
    "deterministic": ["corr", "rmse", "bias"],
}

FAMILY_INPUTS = {
    "prob_binary": ["hcst_prob", "obs_cat", "clim_prob", "n_members"],
    "categorical": ["hcst_prob", "obs_cat", "clim_prob", "n_members"],
    "ensemble_obs": ["obs", "hcst"],
    "deterministic": ["obs", "hcst"],
}

# score: (long_name, units, description). long_name and units label plots; DATA = the units
# of the variable (from the obs' 'units' attribute)
DATA = "data"
DESCRIPTIONS = {
    "bs": ("Brier score", "-", "fair Brier score of the event (lower is better)"),
    "bss": ("Brier skill score", "-", "Brier skill score vs climatological probability"),
    "roc": ("ROC area", "-", "ROC area of the event (= 2AFC for a binary event; 0.5 = no skill)"),
    "rps": ("Ranked probability score", "-", "fair ranked probability score (lower is better)"),
    "rpss": ("Ranked probability skill score", "-",
             "RPS skill score vs climatological probabilities"),
    "heidke": ("Heidke skill score", "-",
               "Heidke skill score of the most likely category (0 = no skill)"),
    "2afc": ("2AFC score", "-", "2AFC score of the category probabilities (0.5 = no skill)"),
    "crps": ("CRPS", DATA, "fair continuous ranked probability score (lower is better)"),
    "crpss": ("CRPS skill score", "-", "CRPS skill score vs leave-one-year-out climatology"),
    "groc": ("Generalized discrimination score", "-",
             "generalized discrimination score of the ensemble (0.5 = no skill)"),
    "corr": ("Correlation", "-", "Pearson correlation, ensemble mean vs obs"),
    "rmse": ("RMSE", DATA, "root-mean-square error of the ensemble mean"),
    "bias": ("Bias", DATA, "mean(ensemble mean - obs)"),
}


# =========================================================================== preparation
def _drop_member_copies(obs):
    return obs.isel(member=0, drop=True) if "member" in obs.dims else obs


def make_categories(obs, hcst, quantiles=None, *, thresholds=None, edges_from="obs"):
    """
    Inputs for the categorical family (and, per category, for prob_binary): categories
    bounded by climatological quantiles, e.g. [1/3, 2/3] (terciles, the default) or
    [0.1, 0.9] (below 10th / normal / above 90th percentile), or by absolute thresholds,
    e.g. thresholds=[1.0] (in the units of the data, e.g. dry / wet weeks for weekly totals).
    edges_from: "obs" (bias-corrected data) or "model" (raw data); see categories.Categories.

    Returns dict(hcst_prob, obs_cat, clim_prob, n_members); all but n_members have the
    dim "category":
        hcst_prob : fraction of (non-NaN) members in each category
        obs_cat   : 1 in the observed category, 0 elsewhere (NaN where obs missing)
        clim_prob : climatological probability of each category (quantile spacing, or
                    observed frequency for absolute thresholds)
        n_members : number of non-NaN members, for the fair scores
    """
    cats = Categories(quantiles, edges_from, thresholds=thresholds).fit(obs, hcst)
    return dict(hcst_prob=cats.apply(hcst), obs_cat=cats.categorise_obs(obs),
                clim_prob=cats.clim_prob, n_members=hcst.count("member"))


def make_event(obs, hcst, quantile=0.9, side="above", *, threshold=None, edges_from="obs"):
    """
    Inputs for the prob_binary family: the event "above the q-th percentile" (side="above")
    or "at or below it" (side="below"); with `threshold` (in the units of the data) the event is
    above / at or below that value instead, and `quantile` is ignored. Same as
    make_categories with one threshold, keeping one category; edges_from as there.
    Returns dict(hcst_prob, obs_cat, clim_prob, n_members) without the category dim.
    """
    if side not in ("above", "below"):
        raise ValueError("side must be 'above' or 'below'")
    k = 1 if side == "above" else 0
    edges = dict(quantiles=[quantile]) if threshold is None else dict(thresholds=[threshold])
    inputs = make_categories(obs, hcst, edges_from=edges_from, **edges)
    return {name: (v.isel(category=k, drop=True) if "category" in v.dims else v)
            for name, v in inputs.items()}


# =========================================================================== helpers
def crps_fair(obs, fc):
    """Fair ensemble CRPS per case; NaN members are ignored."""
    m = fc.count("member")
    term1 = abs(fc - obs).mean("member")
    pair = abs(fc - fc.rename(member="_m2")).sum(["member", "_m2"])
    return term1 - pair / (2 * m * (m - 1))


def loo_climatology_ensemble(obs):
    """For each hindcast year (valid_time), the obs of all other years as an ensemble."""
    dim = "valid_time"
    n = obs.sizes[dim]
    parts = []
    for i in range(n):
        ens = obs.drop_isel({dim: i}).rename({dim: "member"})
        ens = ens.assign_coords(member=np.arange(n - 1))
        parts.append(ens.drop_vars([c for c in ens.coords if dim in ens[c].dims], errors="ignore"))
    return xr.concat(parts, dim=obs[dim]).transpose(dim, "member", ...)


def groc_score(obs, hcst):
    """
    Generalized discrimination score (Weigel & Mason 2011, MWR 139; Mason & Weigel 2009,
    MWR 137). For every pair of years with different observations: F = P(member of one year >
    member of the other) over all member pairs, ties 1/2. Correct if F > 0.5 in the direction
    of the observations (1/2 credit if F = 0.5). 0.5 = no skill, 1 = perfect. Pairs with tied
    observations (e.g. two dry periods) are not counted. One member: (1 + Kendall's tau) / 2.
    """
    dim = "valid_time"
    o = obs.transpose(dim, ...)
    f = hcst.transpose(dim, "member", *o.dims[1:])
    ov = o.values.reshape(o.sizes[dim], -1)                         # (n, P)
    fv = f.values.reshape(f.sizes[dim], f.sizes["member"], -1)    # (n, m, P)
    n = ov.shape[0]
    num = np.zeros(ov.shape[1])
    den = np.zeros(ov.shape[1])
    for i in range(n - 1):
        xi = fv[i]
        for j in range(i + 1, n):
            xj = fv[j]
            gt = xj[:, None, :] > xi[None, :, :]
            eq = xj[:, None, :] == xi[None, :, :]
            valid = np.isfinite(xj)[:, None, :] & np.isfinite(xi)[None, :, :]
            cnt = valid.sum((0, 1))
            F = np.where(cnt > 0, ((gt + 0.5 * eq) * valid).sum((0, 1)) / np.maximum(cnt, 1), np.nan)
            d_obs = np.sign(ov[j] - ov[i])
            use = (d_obs != 0) & np.isfinite(d_obs) & np.isfinite(F)
            credit = np.where(F == 0.5, 0.5, np.where((F > 0.5) == (d_obs > 0), 1.0, 0.0))
            num += np.where(use, credit, 0.0)
            den += use
    score = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    rest = [d for d in o.dims if d != dim]
    return xr.DataArray(score.reshape([o.sizes[d] for d in rest]), dims=rest,
                        coords={d: o[d] for d in rest if d in o.coords})


def two_afc_categorical(obs_cat, hcst_prob):
    """
    2AFC score for probabilistic forecasts of ordered categories (Mason & Weigel 2009,
    MWR 137, section 3d). For every pair of years in different observed categories:
    F = P(category drawn from the forecast of the year observed higher > category drawn
    from the other forecast), not counting draws in the same category. Correct if F > 0.5,
    1/2 credit if F = 0.5 (or if both forecasts are certain of the same category).
    0.5 = no skill, 1 = perfect. With two categories it equals the ROC area.

    obs_cat, hcst_prob : (valid_time, ..., category); obs_cat one-hot (NaN where missing)
    """
    dim = "valid_time"
    o = obs_cat.transpose(dim, "category", ...)
    p = hcst_prob.transpose(*o.dims)
    rest = list(o.dims[2:])
    shape = [o.sizes[d] for d in rest]
    ov = o.values.reshape(o.sizes[dim], o.sizes["category"], -1)            # (n, K, P)
    pv = p.values.reshape(ov.shape)
    valid = np.isfinite(ov).all(1) & np.isfinite(pv).all(1)                  # (n, P)
    cat = np.where(valid, np.nan_to_num(ov).argmax(1), -1)                   # observed category
    cdf_below = np.cumsum(pv, axis=1) - pv                                   # P(X < k), (n, K, P)

    num = np.zeros(ov.shape[2])
    den = np.zeros(ov.shape[2])
    for i in range(ov.shape[0] - 1):
        for j in range(i + 1, ov.shape[0]):
            use = valid[i] & valid[j] & (cat[i] != cat[j])
            hi, lo = np.where(cat[j] > cat[i], j, i), np.where(cat[j] > cat[i], i, j)
            ph = np.where(hi == j, pv[j], pv[i])                             # forecast, higher obs
            cl = np.where(lo == i, cdf_below[i], cdf_below[j])               # other forecast's cdf
            pl = np.where(lo == i, pv[i], pv[j])
            greater = (ph * cl).sum(0)                                       # P(X_hi > X_lo)
            equal = (ph * pl).sum(0)                                         # P(X_hi = X_lo)
            with np.errstate(invalid="ignore", divide="ignore"):
                F = np.where(equal < 1, greater / (1 - equal), 0.5)
            credit = np.where(np.isclose(F, 0.5), 0.5, (F > 0.5).astype(float))
            num += np.where(use, credit, 0.0)
            den += use
    score = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    return xr.DataArray(score.reshape(shape), dims=rest,
                        coords={d: o[d] for d in rest if d in o.coords})


def _family(score):
    for family, names in SKILL_SCORES.items():
        if score in names:
            return family
    raise ValueError(f"unknown score '{score}'; available: {SKILL_SCORES}")


# =========================================================================== main
@quiet_nan_warnings()
def hindcast_skill(score, **inputs):
    """
    Skill of a hindcast, per lead time and grid point.

    score  : name (see SKILL_SCORES) -> DataArray; or list of names -> Dataset
    inputs : the inputs of the score's family (see module docstring), e.g. **make_event(...)
    Scores are computed over valid_time (the hindcast years).
    """
    dim = "valid_time"
    if not isinstance(score, str):
        return xr.Dataset({s: hindcast_skill(s, **inputs)
                           for s in score})

    family = _family(score)
    missing = [k for k in FAMILY_INPUTS[family] if k not in inputs]
    if missing:
        raise ValueError(f"'{score}' ({family}) needs inputs {FAMILY_INPUTS[family]}; "
                         f"missing {missing}")

    # coordinates along the scored dim (e.g. target_time) mean nothing after scoring, and
    # some xskillscore functions carry them through in a garbled form: drop them first
    inputs = {k: v.drop_vars([c for c in v.coords if c not in v.dims and dim in v[c].dims])
              if isinstance(v, xr.DataArray) else v
              for k, v in inputs.items()}

    # ---------------------------------------------------------------- binary event
    if family == "prob_binary":
        p, o = inputs["hcst_prob"], inputs["obs_cat"]
        c, m = inputs["clim_prob"], inputs["n_members"]
        if score in ("bs", "bss"):
            # fair Brier: subtract the expected error due to the finite ensemble
            fair = (p * (1 - p) / (m - 1)).where(m > 1).where(o.notnull())
            bs = xs.brier_score(o, p, dim=dim) - fair.mean(dim)
            if score == "bs":
                out = bs
            else:
                ref = xs.brier_score(o, (xr.zeros_like(p) + c).where(p.notnull()), dim=dim)
                out = (1 - bs / ref).where(ref > 0)        # ref 0: event never/always observed
        elif score == "roc":
            # ROC is undefined unless the event both occurred and did not occur
            # (xs.roc returns arbitrary values there); NaN for missing data as well
            valid = (o.notnull().all(dim) & p.notnull().all(dim)
                     & (o.sum(dim) > 0) & ((1 - o).sum(dim) > 0))
            # xs.roc needs boolean obs without NaN: fill, then mask cells with missing data
            out = xs.roc(o.fillna(0).astype(bool), p.fillna(0), bin_edges="continuous",
                         dim=dim, return_results="area").where(valid)

    # ---------------------------------------------------------------- categories
    elif family == "categorical":
        p, o = inputs["hcst_prob"], inputs["obs_cat"]
        c, m = inputs["clim_prob"], inputs["n_members"]
        if score in ("rps", "rpss"):
            Fc = p.cumsum("category")
            fair = (Fc * (1 - Fc)).sum("category") / (m - 1)
            fair = fair.where(m > 1).where(o.notnull().all("category"))
            rps = xs.rps(o, p, category_edges=None, input_distributions="p", dim=dim) - fair.mean(dim)
            if score == "rps":
                out = rps
            else:
                ref = xs.rps(o, (xr.zeros_like(p) + c).where(p.notnull()), category_edges=None,
                             input_distributions="p", dim=dim)
                out = (1 - rps / ref).where(ref > 0)
        elif score == "heidke":
            # forecast category = most likely one (ties -> the lower category)
            valid = o.notnull().all("category") & p.notnull().all("category")
            fc_cat = p.fillna(-1).argmax("category").where(valid)
            ob_cat = o.fillna(-1).argmax("category").where(valid)
            edges = np.arange(p.sizes["category"] + 1) - 0.5
            out = xs.Contingency(ob_cat, fc_cat, edges, edges, dim=dim).heidke_score()
            out = out.where(valid.any(dim))
        elif score == "2afc":
            out = two_afc_categorical(o, p)

    # ---------------------------------------------------------------- full ensemble
    elif family == "ensemble_obs":
        obs = _drop_member_copies(inputs["obs"])
        hcst = inputs["hcst"]
        if score in ("crps", "crpss"):
            crps = crps_fair(obs, hcst).mean(dim)
            if score == "crps":
                out = crps
            else:
                ref = crps_fair(obs, loo_climatology_ensemble(obs)).mean(dim)
                out = (1 - crps / ref).where(ref > 0)      # ref 0: obs the same every year
        elif score == "groc":
            out = groc_score(obs, hcst)

    # ---------------------------------------------------------------- ensemble mean
    elif family == "deterministic":
        obs = _drop_member_copies(inputs["obs"])
        em = inputs["hcst"].mean("member")
        if score == "corr":
            out = xs.pearson_r(obs, em, dim=dim)
        elif score == "rmse":
            out = xs.rmse(obs, em, dim=dim)
        elif score == "bias":
            out = (em - obs).mean(dim)

    # keep only coordinates that still make sense (e.g. category_label, not target_time)
    out = out.drop_vars([v for v in out.coords if not set(out[v].dims) <= set(out.dims)
                         or (v not in out.dims and not out[v].dims)])
    # attributes from the table only (inherited ones, e.g. the obs' long_name, would
    # mislabel plots)
    long_name, units, text = DESCRIPTIONS[score]
    if units == DATA:
        units = inputs["obs"].attrs.get("units", "")
    out = order_dims(out).rename(score)
    out.attrs = dict(long_name=long_name, units=units, description=text, family=family)
    return out
