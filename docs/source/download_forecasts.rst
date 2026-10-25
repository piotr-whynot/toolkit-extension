Downloading Operational Forecasts
=================================

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.

Download a forecast
-------------------

To download a forecast from ECMWF's S2S database, you will need to use the `download_forecast` function from the `download_forecast.py` python module. To begin with, it is recommended that you download ECMWF forecast data first to get used to the system. 

To import the necessary function, use the following line of python code:

.. code-block:: python
    
    from acacia_s2s_toolkit.download_forecast import download_forecast

After this, use `download_forecast` to download operational forecasts:

.. code-block:: python

   download_forecast(variable, model='ECMWF', fcdate=None, # key forecast parameters
                     start_lt=0, end_lt=None, period='1D', # variables that define leadtimes to download
                     country_name=None, region_name=None, bbox_bounds=[90, -180, -90, 180], # variables that enable you to define spatial grid
                     plevs=None, # optional variable when downloading pressure-level data
                     filename=None, data_save_dir=None, data_format="netcdf", grid="1.5/1.5", # characteristics associated with saved file
                     fc_enslags=None, # defined lags for lagged ensemble
                     overwrite=False, verbose=True, cleanup=True) # options associated with saved outputs 

:Necessary parameters:

- **variable** (*str*): The forecasted variable. Please use variable abbreviations listed on `ECMWF's S2S parameter page <https://confluence.ecmwf.int/display/S2S/Parameters>`_. A few examples include:

  - ``10u``: 10 metre u-velocity (:math:`\mathrm{m \, s^{-1}}`)
  - ``10v``: 10 metre v-velocity (:math:`\mathrm{m \, s^{-1}}`)
  - ``2t``: Surface air temperature (K)
  - ``tp``: Total precipitation (mm)

.. note::

   Each ``variable`` must be requested separately.

:Optional parameters:

- **model** (*str*): The forecasting model. If not specified, defaults to ``'ECMWF'``. A full summary of models including their forecast frequency and data access delay can be found on `ECMWF's Model page <https://confluence.ecmwf.int/display/S2S/Models>`_.  Supported models (as of 6th November 2025) include:

  - ``'ECMWF'``: European Centre for Medium-Range Weather Forecasts (default)
  - ``'ECCC'``: Environment and Climate Change Canada
  - ``'HMCR'``: Hydrometeorological Centre of Russia
  - ``'JMA'``: Japan Meteorological Agency
  - ``'KMA'``: Korea Meteorological Administration
  - ``'NCEP'``: National Centers for Environmental Prediction (NOAA/USA)
  - ``'CMA'``: China Meteorological Administration

- **fcdate** (*str*, optional): Forecast initialisation date in ``YYYYMMDD`` format. If no date is given, then the latest avaliable forecast is downloaded.

- **start_lt** (*int or str*, optional): First lead time (in hours) to include in the downloaded forecast. For example, ``24`` starts from day 1, while ``168`` starts from week 2. Alternatively, a day of the week (e.g. ``"monday"``, ``"tuesday"``) can be provided. In this case, the lead time is automatically calculated so that the forecast starts on the next occurrence of the requested day relative to ``fcdate``. If not specified, defaults to ``0``.

- **end_lt** (*int*, optional): Final lead time (in hours) to include in the downloaded forecast. For example, ``672`` corresponds to day 28 of the forecast. If not requested, set to end of forecast. 

- **period** (*str*, optional): Aggregation period used to select forecast lead times and compute aggregated products. The period should be specified in days, for example ``'1D'``, ``'2D'``, ``'7D'`` or ``'14D'``. Accumulated variables retain only the end of each aggregation period, while instantaneous and averaged variables are aggregated over the specified period. If not requested, set to ``'1D'``.

.. note::

    If no leadtime information provided, the download will focus on **daily** aggregation for the **full forecast length** i.e. ``start_lt=0``, ``end_lt=None``, ``period='1D'``

- **country_name** (*str*, optional): Country name used to derive bounds automatically from Natural Earth country polygons. Overrides ``bbox_bounds`` and ``region_name`` if provided. Current supported values include:

  - ``"ethiopia"``
  - ``"kenya"``
  - ``"madagascar"``
  - ``"eastafrica"``

- **region_name** (*str*, optional): Predefined domain name. Overrides ``bbox_bounds`` if provided. Current supported values include:

  - ``"westafrica"``
  - ``"eastafrica"``
  - ``"southernafrica"``

- **bbox_bounds** (*list of float*, optional, default ``[90, -180, -90, 180]``): Geographic bounding box in the order ``[north, west, south, east]`` (degrees latitude/longitude). Overridden if ``location_name`` is provided.

- **plevs** (*int or list of int*, optional): Pressure levels in hPa for pressure-level variables. Avaliable levels include 1000, 925, 850, 700, 500, 300, 200, 100, 50 and 10 hPa.  

- **filename** (*str*, optional): Name of the outputted file. If ``None``, filename is automatically constructed from variable, model, date, pressure levels, and domain. The function appends ``.nc`` or ``.grib`` depending on ``data_format``.

- **data_save_dir** (*str*, optional): Directory to save the downloaded file. Created if it doesn’t exist. Defaults to the current working directory if ``None``.

- **data_format** (*str*, optional, default ``"netcdf"``): Output data format: ``"netcdf"`` or ``"grib"``. Current testing has focused on netCDF format.

- **grid** (*str*, optional, default ``"1.5/1.5"``): Horizontal grid spacing (degrees) in ``lat/lon`` order.

- **fc_enslags** (*list of int*, optional): The selection of lagged ensemble members relative to the forecast initialisation date. Default values for each model are described in the `ECMWF variables page <https://confluence.ecmwf.int/display/~ecm0847/acacia_s2s_toolkit+available+forecasting+systems>`_.

- **overwrite** (*bool*, default ``False``): If ``True``, overwrite existing file. If ``False``, skip download if file exists.

- **verbose** (*bool*, default ``True``): If ``True``, prints download info, debug messages, and rollback messages. If ``False``, suppresses output.

:Returns:
  Path to the downloaded file as a string.

.. note::

    For simplicity, all ensemble members are downloaded in a single request. Currently, the control and perturbed forecasts are concatenated into one file after download.

You can check the status of your `webAPI downloads here <https://apps.ecmwf.int/webmars/joblist/>`_ . 




