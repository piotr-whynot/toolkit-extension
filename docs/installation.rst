Installation
============

Requirements
------------

toolkit-extension requires the following libraries to be installed manually
via conda before installing the package itself. These libraries have compiled
dependencies that cannot be managed by pip.

.. code-block:: bash

    mamba install -c conda-forge xarray xesmf esmpy netcdf4 h5netcdf

Installing the Package
----------------------

Once the dependencies are installed, install toolkit-extension via pip:

.. code-block:: bash

    pip install toolkit-extension

Verifying the Installation
--------------------------

To verify the installation was successful:

.. code-block:: python

    import toolkit_extension
    print(toolkit_extension.__version__)

Development Installation
------------------------

If you want to contribute or modify the code, clone the repository and install
in editable mode:

.. code-block:: bash

    git clone https://github.com/you/toolkit-extension.git
    cd toolkit-extension
    pip install -e .

Tested Environments
-------------------

toolkit-extension has been tested with the following:

- Python 3.11
- xarray 2024.x
- xesmf 0.8.x
