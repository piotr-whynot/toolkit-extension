rf_enslags (reforecast lagged ensemble)
==================================================

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.

As described in `ECMWF's S2S forecast models <https://confluence.ecmwf.int/display/S2S/Models>`_, different forecasting centres run integrations at varying frequencies and ensemble sizes. Consequently, certain choices must be made to ensure a sufficient number of reforecasts are downloaded.

The ``rf_enslags`` option allows the user to download a **reforecast lagged ensemble**.

- The list of integers provided to create a **reforecast** ensemble can be **less than, equal to, or more than zero**.
- Creating a lagged ensemble will only work if a forecast was initialised on the requested lagged date.

.. note::
   When building a lagged ensemble, the toolkit automatically trims the start and end portions that are not common to every forecast, leaving only the overlapping period where all ensemble members are available.  This is described in more detail in the `Lagged Ensembles <https://acacia-s2s-toolkit.readthedocs.io/en/latest/theory_lagged_ensembles.html>`_ theory.
