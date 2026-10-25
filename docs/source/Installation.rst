Getting started
===============

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.

.. _installation:

Installation
------------

To use the **acacia-s2s-toolkit**, you will first need to have `Python 3 <https://www.python.org/downloads/>`_ installed.

You can then install the package using ``pip``:

.. code-block:: console

   pip install acacia_s2s_toolkit

Upgrading the Package
---------------------
To upgrade to the latest version, run:

.. code-block:: bash

   python3 -m pip install --upgrade acacia_s2s_toolkit

.. note::
    This tool is being developed by the ACACIA project. Regular package updates are recommended due to continuous developments. 

Necessary requirements
----------------------

Saving ECMWF Data Store credentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
At the time of writing, sub-seasonal forecast data from `ECMWF's S2S database <https://ecds.ecmwf.int/datasets/s2s-forecasts?tab=overview>`_ is made available through ECMWF's Data Store. To download data from the ECMWF Data Store, you need to configure your CDS API key.

After registering for an ECMWF account, your API credentials can be found `here <https://ecds.ecmwf.int/how-to-api>`_. 

For the acacia_s2s_toolkit, please save these credentials to your home directory as a file named ``~/.cdsapirc.ecds``. 

.. important::
    Please note that the credentials should be saved in a file called ``~/.cdsapirc.ecds`` and not the default ``~/.cdsapirc``. This is to avoid confusion with API credentials accessing ECMWF Climate Data Store (i.e. for ERA5 data). 

Installing ECMWF client libraries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Once your ECMWF webAPI credentials are saved, install the CDS API client using ``pip``:

.. code-block:: console

   pip install "cdsapi>=0.7.7"

This package provides the necessary client functions for accessing ECMWF data.

Additional package dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In addition to requiring "*cdsapi*", this package depends on the following libraries:

- **numpy** (version 1.23 or higher)
- **cdo** (versions 2.4.0 or higher)
- **xarray** (version 2024.09.0 or higher)
- **eccodes** (version 2.40.0 or higher)
- **dask** (version 2024.9.0)
- **pandas** (version 2.2.3 or higher)
- **scipy** (version 1.14.1 or higher)
- **netCDF4** (version 1.7.2 or higher)
- **requests** (versions 2.32.2 or higher)
- **matplotlib** (versions 3.8 or higher)
- **cartopy** (versions 0.22 or higher)

.. note::
    Some dependencies, especially eccodes, cdo, and cartopy—require system-level libraries (e.g., PROJ, GEOS, ecCodes) in addition to the Python packages. Using Conda is recommended for effortless installation of both Python and non-Python components.

Given the number of dependencies, we highly recommend using the **acacia_s2s_toolkit** with a dedicated virtual environment. Guidance on creating and managing environments can be found in the official documentation for both `Conda <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ and `Anaconda <https://www.anaconda.com/docs/getting-started/working-with-conda/environments>`_. 

Quick checklist
~~~~~~~~~~~~~~~
Before attempting to download sub-seasonal forecast data using this python tool, please ensure you have performed the following steps:

1. `Python 3 <https://www.python.org/downloads/>`_ is installed
2. The Python packages **acacia_s2s_toolkit**, **cdsapi** and **sites-toolkit** are installed. 
3. Your ECMWF webAPI credentials are saved in ``~/.cdsapirc.ecds`` on your working environment.


Further information
-------------------

For more details on setting up and using ECMWF Data Store, please visit `using the ECMWF data store <https://confluence.ecmwf.int/display/DAC/Using+the+ECMWF+data+store+-+ECDS>`_

