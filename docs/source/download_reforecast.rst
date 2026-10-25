Downloading Operational Reforecasts
===================================

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.

Download a collection of reforecasts
------------------------------------

To download a selection of reforecasts from ECMWF's S2S database, you will need to use the `download_hindcast` function from the `download_hindcast.py` python module. To begin with, it is recommended that you download ECMWF reforecast data first to get used to the system. 

.. note::

   In this toolkit, the words 'reforecast' and 'hindcast' can be used interchangeably.  The function names tend to use the convention 'hindcast', whilst the documentation uses the convention 'reforecast'.

To import the necessary function, use the following line of python code:

.. code-block:: python
    
    from acacia_s2s_toolkit.download_hindcast import download_hindcast

After this, use the `download_hindcast` function to download operational reforecasts:

.. code-block:: python

   download_hindcast(variable, model='ECMWF', fcdate=None, # key forecast parameters
                     start_lt=0, end_lt=None, period='1D', # variables that define leadtimes to download
                     country_name=None, region_name=None, bbox_bounds=[90, -180, -90, 180], # variables that enable you to define spatial grid
                     plevs=None, # optional variable when downloading pressure-level data
                     filename=None, data_save_dir=None, data_format="netcdf", grid="1.5/1.5", # characteristics associated with saved file
                     rf_enslags=None, rf_years=None, fc_time=True, # defined lags for lagged ensemble + time coordinate
                     overwrite=False, verbose=True, cleanup=True) # options associated with saved outputs

This function retrieves operational reforecast (hindcast) data for a specified model and variable from the S2S database. It supports downloading multiple initialization years and ensemble lags for a given forecasting system.

:Necessary parameters:

- **variable** (*str*): The forecasted variable to download. Use variable abbreviations listed on `ECMWF's S2S parameter page <https://confluence.ecmwf.int/display/S2S/Parameters>`_.

.. note::

   Each ``variable`` must be requested separately.

:Optional parameters:

Parameters including **model**, **fcdate**, **start_lt**, **end_lt**, **period**, **country_name**, **region_name**, **bbox_bounds**, **plevs**, **filename**, **data_save_dir**, **data_format**, **grid**, **overwrite** and **verbose** are described in detail on the `download_forecasts <https://acacia-s2s-toolkit.readthedocs.io/en/latest/download_forecasts.html>`_ webpage. 

.. note::

   The reforecast should use the same ``model`` as the forecast.

Additional parameters specific to reforecast downloads include:

- **rf_years** (*list of int*, optional): List of years to download reforecasts for. For example, [2000, 2001, 2002] will retrieve all reforecasts initialised between 2000 and 2002. Please note, due to changes when transferring to ECMWF Climate Data Store, all rf_years will be downloaded between the minimum and maximum of requested rf_years. 

- **rf_enslags** (*list of int*, optional): List of lagged ensemble members (in days) relative to the forecast initialization date. Default values depend on the selected model; see the following `confluence page <https://confluence.ecmwf.int/display/~ecm0847/acacia_s2s_toolkit+available+forecasting+systems>`_ for details. 

- **fc_time** (*bool*, default ``True``): **Option in development** If ``True``, shifts all lagged ensemble members so their forecast times are aligned to a common initialization date (useful for ensemble averaging). If False, preserves each member’s original valid time (useful for forecast validation).

:Returns:
  Path to the downloaded file as a string.

.. note::

    For simplicity, all ensemble members are downloaded in a single request. Currently, the control and perturbed reforecasts are concatenated into one file after download.

You can check the status of your webAPI downloads on the following `page <https://apps.ecmwf.int/webmars/joblist/>`_. 



