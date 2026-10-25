"""
Categories defined by thresholds: fit once, apply to any ensemble or to obs.

    cats = Categories([1/3, 2/3]).fit(obs, hcst)         # quantiles (default: terciles)
    cats = Categories(thresholds=[1.0]).fit(obs, hcst)   # absolute values, in data units
    hc_prob = cats.apply(hcst)                           # (valid_time, effective_lead_time, category, latitude, longitude)
    fc_prob = cats.apply(fcst)                           # same thresholds
    obs_cat = cats.categorise_obs(obs)                   # 1 in the observed category, else 0

    cats.to_dataset() / Categories.from_dataset(ds)      # thresholds saved with the product

Category k contains values with threshold[k-1] < x <= threshold[k]. With quantiles [0.9]
there are two categories: 0 = at or below the 90th percentile, 1 = above it; with
thresholds [1.0]: 0 = at most 1, 1 = above 1.

Quantile thresholds are computed over all hindcast years, per lead time and grid point:
    edges_from="obs"   : obs quantiles for both obs and ensembles (bias-corrected data)
    edges_from="model" : obs quantiles for obs, quantiles of the hindcast (years and members
                         pooled) for ensembles (raw model data)
Absolute thresholds are the same everywhere and for obs and ensembles; use bias-corrected
data, since a model bias shifts the forecast probabilities.

Climatological probability of each category (the reference forecast for skill scores):
    quantiles : the quantile spacing, e.g. 1/3 each for terciles
    absolute  : the observed frequency over the hindcast years, per lead time and grid point
"""

from __future__ import annotations

import warnings
from contextlib import contextmanager

import numpy as np
import xarray as xr

TERCILES = [1 / 3, 2 / 3]


# numpy / xskillscore warnings caused by cells without data (e.g. ocean), where the result
# is NaN anyway; only these are silenced
_NAN_WARNINGS = ("All-NaN slice", "Mean of empty slice", "Degrees of freedom",
                 "invalid value encountered", "divide by zero encountered",
                 "Converting input from bool")


@contextmanager
def quiet_nan_warnings():
    """Context manager / decorator silencing the warnings in _NAN_WARNINGS."""
    with warnings.catch_warnings(), np.errstate(invalid="ignore", divide="ignore"):
        for msg in _NAN_WARNINGS:
            warnings.filterwarnings("ignore", message=msg, category=RuntimeWarning)
        yield


def order_dims(da):
    """latitude, longitude last (CF order Y, X); other dims, e.g. category, just before them,
    like member: (valid_time, effective_lead_time, category, latitude, longitude)."""
    return da.transpose(..., "latitude", "longitude", missing_dims="ignore")


def _drop_member_copies(obs):
    return obs.isel(member=0, drop=True) if "member" in obs.dims else obs


def _labels(edges, prefix):
    e = [f"{prefix}{v:.3g}" for v in edges]
    return [f"<={e[0]}"] + [f"{a}-{b}" for a, b in zip(e[:-1], e[1:])] + [f">{e[-1]}"]


class Categories:
    """Quantile-based or absolute-threshold categories; see module docstring."""

    def __init__(self, quantiles=None, edges_from="obs", *, thresholds=None):
        if quantiles is not None and thresholds is not None:
            raise ValueError("give quantiles or thresholds, not both")
        if edges_from not in ("obs", "model"):
            raise ValueError("edges_from must be 'obs' or 'model'")
        if thresholds is not None:
            if edges_from == "model":
                raise ValueError("edges_from='model' applies to quantiles only; absolute "
                                 "thresholds are the same for obs and model")
            self.kind, edges = "absolute", thresholds
        else:
            self.kind, edges = "quantile", TERCILES if quantiles is None else quantiles
        e = [float(v) for v in np.atleast_1d(edges)]
        if not all(a < b for a, b in zip(e, e[1:])):
            raise ValueError("thresholds must be increasing")
        if self.kind == "quantile" and not (0 < e[0] and e[-1] < 1):
            raise ValueError("quantiles must be between 0 and 1")
        self.edges, self.edges_from = e, edges_from

    # ------------------------------------------------------------------ fit
    @quiet_nan_warnings()
    def fit(self, obs, hcst):
        """Thresholds (and, for absolute thresholds, observed frequencies) from the hindcast
        period."""
        obs = _drop_member_copies(obs)
        if self.kind == "quantile":
            self.threshold_obs = obs.quantile(self.edges, dim="valid_time").rename(quantile="edge")
            if self.edges_from == "obs":
                self.threshold_model = self.threshold_obs
            else:
                self.threshold_model = hcst.quantile(self.edges, dim=["valid_time", "member"]) \
                                           .rename(quantile="edge")
            p = np.diff(np.r_[0.0, self.edges, 1.0])
            self._clim = xr.DataArray(p, dims="category", coords=self._category_coords())
        else:
            # the same values at every lead time and grid point (keeps the grid check)
            grid = xr.ones_like(obs.isel(valid_time=0, drop=True))
            grid = grid.drop_vars([c for c in grid.coords if c not in grid.dims])
            values = xr.DataArray(self.edges, dims="edge", coords={"edge": self.edges})
            self.threshold_obs = self.threshold_model = values * grid
            self._clim = order_dims(self.categorise_obs(obs).mean("valid_time"))
        return self

    @property
    def clim_prob(self):
        """Climatological probability of each category (see module docstring)."""
        return self._clim

    # ------------------------------------------------------------------ apply
    def apply(self, ens):
        """Fraction of (non-NaN) members in each category; dim 'category' replaces members."""
        return order_dims(self._one_hot(ens, self.threshold_model).mean("member"))

    def categorise_obs(self, obs):
        """1 in the observed category, 0 in the others, NaN where obs are missing."""
        obs = _drop_member_copies(obs)
        return order_dims(self._one_hot(obs, self.threshold_obs))

    # ------------------------------------------------------------------ save / load
    def to_dataset(self):
        """Thresholds and settings as a Dataset (to be stored in the product file)."""
        ds = xr.Dataset({"threshold_obs": self.threshold_obs, "clim_prob": self.clim_prob})
        if self.edges_from == "model":
            ds["threshold_model"] = self.threshold_model
        ds["threshold_obs"].attrs["description"] = "category thresholds for observations"
        ds["clim_prob"].attrs["description"] = "climatological probability of each category"
        return ds.assign_attrs(threshold_type=self.kind, thresholds=self.edges,
                               edges_from=self.edges_from)

    @classmethod
    def from_dataset(cls, ds):
        """Rebuild fitted Categories from to_dataset() output (e.g. a saved product file)."""
        a = ds.attrs
        edges = list(np.atleast_1d(a["thresholds"]))
        kw = dict(quantiles=edges) if a["threshold_type"] == "quantile" else dict(thresholds=edges)
        cats = cls(edges_from=a["edges_from"], **kw)
        cats.threshold_obs = ds["threshold_obs"]
        cats.threshold_model = ds["threshold_model"] if a["edges_from"] == "model" else ds["threshold_obs"]
        cats._clim = ds["clim_prob"]
        return cats

    # ------------------------------------------------------------------ internals
    def _category_coords(self):
        labels = _labels(self.edges, "q" if self.kind == "quantile" else "")
        return {"category": np.arange(len(self.edges) + 1), "category_label": ("category", labels)}

    def _one_hot(self, x, thr):
        if not hasattr(self, "threshold_obs"):
            raise RuntimeError("Categories not fitted; call fit() or from_dataset() first")
        try:   # same lead times and grid as the data the thresholds were fitted on
            xr.align(x, thr, join="exact", exclude=["valid_time", "hcst_valid_time", "member", "edge"])
        except ValueError as err:
            raise ValueError("data and thresholds have different coordinates (lead times or "
                             f"grid); were they prepared the same way?\n{err}") from None
        k = xr.DataArray(np.arange(len(self.edges) + 1), dims="category")
        idx = (x > thr).sum("edge")                       # category index of every value
        out = (idx == k).where(x.notnull())
        return out.assign_coords(self._category_coords())
