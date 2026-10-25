fc_enslags (forecast lagged ensemble)
=======================================

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.

As described in `ECMWF's S2S forecast models <https://confluence.ecmwf.int/display/S2S/Models>`_, different forecasting centres run integrations at varying frequencies and ensemble sizes. Consequently, certain choices must be made to ensure a sufficient number of forecasts are downloaded.

The ``fc_enslags`` option allows the user to download a **forecast lagged ensemble**.

- The list of integers provided to create a **forecast** ensemble must always be **less than or equal to zero**.
- Creating a lagged ensemble will only work if a forecast was initialised on the requested lagged date.

.. note::
    When building a lagged ensemble, the toolkit automatically trims the start and end portions that are not common to every forecast, leaving only the overlapping period where all ensemble members are available.  This is described in more detail in the `Lagged Ensembles <https://acacia-s2s-toolkit.readthedocs.io/en/latest/theory_lagged_ensembles.html>`_ theory.

Example 1: Single forecast initialisation
-----------------------------------------

Suppose you want to request a JMA forecast of accumulated precipitation on **6th November 2025**. You can set ``fc_enslags=[0]``:

.. code-block:: python

   download_forecast('JMA', 'tp', fcdate='20251106', fc_enslags=[0])

This will download all timesteps of forecasted precipitation from the forecast run initialised on 6th November 2025. Since ``fc_enslags`` is set to ``0``, it only includes integrations starting on that date.

Example 2: Lagged forecast ensemble
-----------------------------------

If you set ``fc_enslags=[0, -1, -2]`` (the default for JMA forecasts):

.. code-block:: python

   download_forecast('JMA', 'tp', fcdate='20251106', fc_enslags=[0, -1, -2])

This will download forecasts from **three initialisations**:

- 2025-11-06 (``0``)
- 2025-11-05 (``-1``)
- 2025-11-04 (``-2``)

For the lagged ensemble, the **time dimension corresponds to the forecasted period**. For example, when requesting a lag of ``-2``, the first two days of the forecast are removed so that all forecasts align at the same starting point.

.. note::

   To download all available forecast data from a model, it is recommended to use ``fc_enslags=[0]`` and loop over start dates.
