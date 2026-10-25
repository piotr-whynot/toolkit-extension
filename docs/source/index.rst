Welcome to acacia_s2s_toolkit's documentation!
==============================================

**acacia_s2s_toolkit** is a Python library designed to support sub-seasonal operational forecasting and model evaluation. Specifically, this python wrapper is designed to support ACACIA partners, however it can be used by collaborators. 

The S2S Database is a global archive of sub-seasonal forecasts and re-forecasts established by the joint World Weather Research Programme (WWRP) and World Climate Research Programme (WCRP) Sub-seasonal to Seasonal Prediction project.  It is hosted by the UK’s European Centre for Medium-range Weather Forecasts (ECMWF) and can be accessed via the `ECMWF Data Store <https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=overview>`_.  

The python wrapper pulls data from the ECMWF Data Store through creating appropriate request scripts to retrieve data efficiently.

.. note::

   This project is under active development.

.. important::

   This python wrapper is originally designed with operational forecasters in mind. Please be careful with using the python wrapper for forecast evaluation as certain timestamps may be invalid for your purpose.

Main authors
------------

- Joshua Talib (ECMWF), joshua.talib@ecmwf.int
- Innocent Masukwedza (NCAS-Reading), g.t.masukwedza@reading.ac.uk
- Piotr Wolski (University of Cape Town), wolski@csag.uct.ac.za
- Linda Hirons (NCAS-Reading), l.c.hirons@reading.ac.uk

Contents
--------
.. toctree::
   :maxdepth: 1
   :caption: Core Toolkit Functions

   Installation
   download_forecasts
   download_reforecast
   bias_correction
   skill_evaluation
   
.. toctree::
   :maxdepth: 1
   :caption: Understanding User Inputs

   models 
   fc_enslags
   rf_enslags

.. toctree::
   :maxdepth: 1
   :caption: Understanding Toolkit Outputs

   output
   grid

.. toctree::
   :maxdepth: 1
   :caption: Example notebooks

   notebooks/deterministic_forecast_example
   notebooks/probabilistic_forecast_example
   notebooks/bias_correction_example
   notebooks/postprocessing_basic_example

.. toctree::
   :maxdepth: 1
   :caption: Theory

   theory_subseasonal_forecasting
   theory_reforecasts
   theory_lagged_ensembles
   theory_bias_correction
   tips_faq
   

   
   
   
   
