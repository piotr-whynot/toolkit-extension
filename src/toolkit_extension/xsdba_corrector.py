"""
Bias correctors working on data in corrector format (made by bias_correction._to_bc):

    ref  : (time, effective_lead_time, ...)           observations, no member dimension
    hist : (time, effective_lead_time, member, ...)   hindcast; no member-indexed coordinates
    sim  : (time, effective_lead_time, member, ...)   hindcast or forecast; any number of members

'time' holds the valid_time values (xsdba requires this name); the remaining dimensions
(latitude/longitude or point) are kept separate: each gets its own transfer function.
"""

from __future__ import annotations

import copy
import inspect
import warnings
from abc import ABC, abstractmethod

import numpy as np
import xarray as xr
import xsdba


class BiasCorrector(ABC):
    """
    Common interface for all bias-correction methods.

    Subclasses must implement fit() and apply() on data in corrector format (module
    docstring). They inherit fit_apply_loo().
    """

    @abstractmethod
    def fit(self, ref: xr.DataArray, hist: xr.DataArray) -> "BiasCorrector": ...

    @abstractmethod
    def apply(self, sim: xr.DataArray) -> xr.DataArray: ...

    def fit_apply_loo(self, ref: xr.DataArray, hist: xr.DataArray, leave_out: int = 1) -> xr.DataArray:
        """
        Cross-validated hindcast: for each group of `leave_out` consecutive years, fit a copy
        on the remaining years and correct the held-out ones. Returns hist-shaped output.
        Works for any subclass, since it only uses fit() and apply().
        Does not change the fitted state of this object.
        """
        years = hist["time"].dt.year.values
        uniq = np.unique(years)
        pieces = []
        for i in range(0, len(uniq), leave_out):
            test = np.isin(years, uniq[i:i + leave_out])
            model = copy.deepcopy(self)          # same settings, its own fitted state
            model.fit(ref.isel(time=~test), hist.isel(time=~test))
            pieces.append(model.apply(hist.isel(time=test)))
        return xr.concat(pieces, dim="time").transpose(*hist.dims)


def _make_lead_blocks(n_lead_times: int, window_size: int) -> list[np.ndarray]:
    """Positional indices of consecutive blocks of `window_size` leads.
    A trailing block shorter than `window_size` is merged into the previous one."""
    if window_size < 1 or n_lead_times < 1:
        raise ValueError("n_lead_times and window_size must be >= 1")
    n_blocks = max(n_lead_times // window_size, 1)
    edges = [i * window_size for i in range(n_blocks)] + [n_lead_times]
    return [np.arange(edges[i], edges[i + 1]) for i in range(n_blocks)]



XSDBA_METHODS = {
    "eqm": xsdba.EmpiricalQuantileMapping,
    "scaling": xsdba.Scaling,
    "dqm": xsdba.DetrendedQuantileMapping,
}

# Options that run but give wrong or unusable results in this setup (see class docstring).
_UNSUPPORTED_PARAMS = {
    "adapt_freq_thresh": "dry-day frequency adaptation is computed per member, not pooled, and "
                         "fails when sim has a different number of members than hist",
    "group": "the Grouper is built by the corrector from `pool_dims` and `pooling`",
}


class XsdbaCorrector(BiasCorrector):
    """
    Bias correction of subseasonal hindcasts/forecasts with xsdba, with three ways of
    pooling data across lead times.

    Works on data in corrector format (module docstring). ``time`` must be datetime and all
    arrays need a ``units`` attribute. float32 input is converted to float64 internally and
    the result is returned in the dtype of ``sim``.

    Parameters
    ----------
    method : {"eqm", "scaling", "dqm"} or xsdba adjustment class, default "eqm"
        Bias-correction method; the options each method accepts are listed under
        *Methods and their options* below. A class can be passed for other univariate
        xsdba methods that accept a ``group`` argument (untested).
    pooling : {"lead", "rolling", "blocks"}, default "lead"
        How data are pooled across lead times. Every transfer function always pools all
        training years (``time``) and ``pool_dims``.

        "lead"     One transfer function per lead time. Use for aggregated data
                   (e.g. weekly means), where each lead is already a block.
        "rolling"  One transfer function per lead time, trained on a centred window of
                   ``window_size`` leads. Corrections change smoothly with lead. Near the
                   first and last leads the window is truncated (fewer samples).
        "blocks"   One transfer function per block of ``window_size`` consecutive leads,
                   applied to every lead in the block. A trailing block shorter than
                   ``window_size`` is merged into the previous one (e.g. 30 leads with
                   window_size=7 -> blocks of 7, 7, 7, 9).
    window_size : int, default 7
        Number of leads per window or block. Ignored for pooling="lead".
    pool_dims : sequence of str, default ("member",)
        Dimensions pooled into each transfer function in addition to ``time``. Add e.g.
        "point" to fit a single transfer function for the whole domain.
    **params
        Options of the chosen method (below). They are routed to ``train`` or ``adjust``
        automatically, using the method's own signatures.

    Methods and their options
    -------------------------
    Sample size per transfer function = n_years x n_members (x window_size for "rolling"
    and "blocks"), but only n_years (x window_size) observations. Choose ``nquantiles``
    with the observation count in mind.

    "eqm" : Empirical Quantile Mapping (xsdba.EmpiricalQuantileMapping)
        Maps each value through the hindcast CDF onto the observed CDF: F_ref^-1(F_hist(x)).
        Corrects mean, spread and shape. Recommended default.

        nquantiles : int or 1-D array of floats in (0, 1), default 20
            Number of quantiles, or the quantile levels themselves. xsdba adds end points
            at 1e-6 and 1 - 1e-6. Keep well below the number of observations per transfer
            function (about 10 for 20 obs).
        kind : {"+", "*"}, default "+"
            Additive (temperature) or multiplicative (precipitation, other positive
            variables) correction.
        interp : {"nearest", "linear", "cubic"}, default "nearest"
            Interpolation of the correction between quantiles. "linear" avoids the step
            function of "nearest". "cubic" can produce negative values for precipitation.
        extrapolation : {"constant", "nan"}, default "constant"
            Values outside the training range get the correction of the outermost
            quantile ("constant") or NaN ("nan"; a forecast often has such values).
        jitter_under_thresh_value : str, optional
            Replace values below this threshold by uniform random noise in (0, threshold)
            before training, e.g. "0.01 mm/d" for precipitation with many zeros (ties
            between quantiles). Training then becomes slightly random between runs.
        max_tail_factor : float, optional
            Values above max_tail_factor x the highest hindcast quantile are left
            uncorrected, guarding against extreme corrections in the tail.

    "scaling" : mean scaling (xsdba.Scaling)
        Adds (kind="+") or multiplies (kind="*") a single factor so that the hindcast mean
        matches the observed mean. Spread and shape are not changed. Robust with few years.

        kind : {"+", "*"}, default "+"
        interp : {"nearest", "linear", "cubic"}, default "nearest"
            Only relevant when xsdba groups within ``time`` (not used here).

    "dqm" : Detrended Quantile Mapping (xsdba.DetrendedQuantileMapping)
        EQM applied to detrended data, with the trend added back afterwards. Designed for
        climate projections. Not recommended for forecasts: with a single forecast date the
        trend of ``sim`` is not defined (xsdba warns "Polyfit may be poorly conditioned").

        nquantiles, kind, interp, extrapolation, jitter_under_thresh_value,
        max_tail_factor : as for "eqm"
        detrend : int or xsdba detrending object, default 1
            Degree of the polynomial trend removed from ``sim``.

    Not supported
    -------------
    adapt_freq_thresh (eqm, dqm)
        Raises ValueError. The dry-day frequency correction is computed per member rather
        than pooled, and adjusting a forecast with a different member count fails.
    group
        Built by the corrector from ``pool_dims`` and ``pooling``.
    QuantileDeltaMapping
        Ranks ``sim`` within its own time series; with one forecast date every value
        becomes NaN.
    LOCI
        Fails when members are pooled (xsdba error "tuple.index(x): x not in tuple").
    ExtremeValues, PrincipalComponents, MBCn, NpdfTransform, OTC, dOTC
        Need inputs beyond (ref, hist, sim) or a multivariate dimension.

    Attributes
    ----------
    models : list of trained xsdba adjustment objects (one, or one per block)
    blocks : list of positional lead indices per block (pooling="blocks" only)
    train_kw, adjust_kw : dict
        The routed options.

    Methods
    -------
    fit(ref, hist)
        Train the transfer function(s). Returns self.
    apply(sim)
        Correct ``sim`` (hindcast or forecast). ``sim`` must have the same
        effective_lead_time values as ``hist``. Returns an array with the shape, dims order
        and dtype of ``sim``.
    fit_apply_loo(ref, hist, leave_out=1)
        Cross-validated hindcast: for each group of ``leave_out`` consecutive years, fit on
        the remaining years and correct the held-out ones. Does not change the fitted state.

    Notes
    -----
    - All-NaN slices (e.g. lagged members that have not started or have ended) give NaN in
      the output; the corresponding numpy warning is silenced.
    - Masked cells (e.g. ocean) are removed beforehand by bias_correct (mask=...).

    Examples
    --------
    Aggregated (weekly) temperature:

    >>> bc = XsdbaCorrector("eqm", pooling="lead", nquantiles=10, kind="+", interp="linear")
    >>> bc.fit(ref, hist)
    >>> fcst_bc = bc.apply(sim_fcst)
    >>> hcst_cv = bc.fit_apply_loo(ref, hist)

    Daily precipitation, weekly blocks:

    >>> bc = XsdbaCorrector("eqm", pooling="blocks", window_size=7, nquantiles=15, kind="*",
    ...                     interp="linear", jitter_under_thresh_value="0.01 mm/d")

    Daily temperature, mean correction in a 7-day rolling window:

    >>> bc = XsdbaCorrector("scaling", pooling="rolling", window_size=7, kind="+")
    """

    POOLINGS = ("lead", "rolling", "blocks")
    _WIN = "lead_win"

    def __init__(self, method="eqm", *, pooling="lead", window_size=7, pool_dims=("member",),
                 **params):
        if pooling not in self.POOLINGS:
            raise ValueError(f"pooling must be one of {self.POOLINGS}")
        if isinstance(method, str) and method not in XSDBA_METHODS:
            raise ValueError(f"method must be one of {list(XSDBA_METHODS)} or an xsdba class")
        bad = set(params) & set(_UNSUPPORTED_PARAMS)
        if bad:
            k = sorted(bad)[0]
            raise ValueError(f"'{k}' is not supported: {_UNSUPPORTED_PARAMS[k]}")
        self.cls = XSDBA_METHODS[method] if isinstance(method, str) else method
        self.pooling, self.window_size = pooling, int(window_size)
        self.pool_dims = list(pool_dims)

        adjust_names = set(inspect.signature(self.cls._adjust).parameters) - {"self", "sim"}
        self.train_kw = {k: v for k, v in params.items() if k not in adjust_names}
        self.adjust_kw = {k: v for k, v in params.items() if k in adjust_names}

        self.models = None   # trained xsdba objects: one, or one per block
        self.blocks = None   # positional lead indices per block
        self._lead = None    # lead coordinate of the training data

    # ------------------------------------------------------------------ checks (no conversion)
    @staticmethod
    def _check(da: xr.DataArray, training: bool, name: str):
        if "time" not in da.dims:
            raise ValueError(f"{name}: no 'time' dimension; convert to the corrector format first")
        if "units" not in da.attrs:
            raise ValueError(f"{name}: no 'units' attribute")
        if training:
            bad = [c for c in da.coords if "member" in da[c].dims]
            if bad:
                raise ValueError(f"{name}: member-indexed coordinates {bad} must be removed "
                                 "for training")

    # ------------------------------------------------------------------ internals
    @staticmethod
    def _f64(da):
        # Compute in float64. With float32 input, scipy's linear interpolation (interp="linear")
        # takes a code path that returns NaN wherever adjacent quantiles are tied
        # (e.g. zeros in precipitation, coarsely resolved data).
        return da if da.dtype == np.float64 else da.astype(np.float64).assign_attrs(da.attrs)

    def _windows(self, da):
        w = da.rolling(effective_lead_time=self.window_size, center=True).construct(self._WIN)
        return w.assign_attrs(da.attrs)

    def _train(self, ref, hist, add_dims):
        group = xsdba.Grouper("time", add_dims=add_dims)
        return self.cls.train(ref, hist, group=group, **self.train_kw)

    def _adjust(self, model, sim, dummy=None):
        # Members that have not started / have finished (lagged ensembles) are all-NaN
        # slices; numpy warns about them, the output is correctly NaN. Silence only that.
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message="All-nan slice encountered",
                                    category=RuntimeWarning)
            if dummy is None:
                return model.adjust(sim, **self.adjust_kw)
            out = model.adjust(sim.expand_dims(dummy), **self.adjust_kw)   # dummy dim, no coordinate
            return out.squeeze(dummy, drop=True)

    # ------------------------------------------------------------------ public API
    def fit(self, ref: xr.DataArray, hist: xr.DataArray) -> "XsdbaCorrector":
        self._check(ref, True, "ref")
        self._check(hist, True, "hist")
        ref, hist = self._f64(ref), self._f64(hist)
        self._lead = hist["effective_lead_time"].values

        if self.pooling == "lead":
            self.models = [self._train(ref, hist, self.pool_dims)]
            self.blocks = None
        elif self.pooling == "rolling":
            self.models = [self._train(self._windows(ref), self._windows(hist),
                                       self.pool_dims + [self._WIN])]
            self.blocks = None
        else:
            self.blocks = _make_lead_blocks(hist.sizes["effective_lead_time"], self.window_size)
            self.models = [self._train(ref.isel(effective_lead_time=b),
                                       hist.isel(effective_lead_time=b),
                                       self.pool_dims + ["effective_lead_time"])
                           for b in self.blocks]
        return self

    def apply(self, sim: xr.DataArray) -> xr.DataArray:
        if self.models is None:
            raise RuntimeError("call fit() first")
        self._check(sim, False, "sim")
        if not np.array_equal(sim["effective_lead_time"].values, self._lead):
            raise ValueError("'effective_lead_time' of sim differs from the training data")
        x = self._f64(sim)

        if self.pooling == "lead":
            out = self._adjust(self.models[0], x)
        elif self.pooling == "rolling":
            out = self._adjust(self.models[0], x, dummy=self._WIN)
        else:
            out = xr.concat([self._adjust(m, x.isel(effective_lead_time=b))
                             for m, b in zip(self.models, self.blocks)],
                            dim="effective_lead_time")

        out = out.transpose(*sim.dims)
        if np.issubdtype(sim.dtype, np.floating) and out.dtype != sim.dtype:
            out = out.astype(sim.dtype).assign_attrs(out.attrs)     # back to the input dtype
        if out.sizes != sim.sizes:
            raise RuntimeError(f"size mismatch after adjust: {dict(out.sizes)} vs {dict(sim.sizes)}")
        return out
