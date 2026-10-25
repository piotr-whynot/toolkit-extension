Bias Correction
===============

.. important::
    This is a tool being developed by the ACACIA project. Please email feedback to joshua.talib@ecmwf.int.
    
For bias correction of forecasts and reforecasts, you will need to use the `bias_correction` function from the `bias_correction.py` python module.  The function uses `xsdba <https://xsdba.readthedocs.io/en/latest/#>`_, a Python library for statistical downscaling and bias adjustment.  

To import the necessary functions, use the following line of python code:

.. code-block:: python
    
    from acacia_s2s_toolkit.bias_correction import bias_correction

After this, use `bias_correction` as follows:

.. code-block:: python

   bias_correction(obs, hcst, fcst, # key bias correction data
		   method="dqm", cv=True, # parameters that define bias correction method
		   mask="auto", units, return_model 
		   options) # additional options for xsdba bias correction methods

:Parameters:

- **obs**: (valid_time, lead, lat, lon) observations are aligned to hcst according to ``align_obs`` function.

- **hcst**: (valid_time, lead, lat, lon)

- **fcst** (optional): (valid_time, lead, lat, lon)

- **method** (*str*): the bias correction method of choice.

  - ``'scaling'``: simple mean bias correction (adjusts model average) using xsdba
  - ``'meanvar'``: correcting both mean and variance (adjusts model average and spread)
  - ``'eqm'``: empirical quantile mapping (all model quantiles adjusted separately) using xsdba
  - ``'dqm'``: detrended quantile mapping (all model quantiles adjusted separately while conserving long-term trends) using xsdba

- **cv**: if ``'True'`` the returned hindcast is cross-validated (leave_out years at a time), as needed for skill assessment.  If ``'False'`` it is corrected in-sample.  The forecast is always corrected with the model fitted on all years.  

- **mask**: "auto" (cells where obs have data), a boolean (lat, lon) array, or None.  Cells outside the mask are not corrected and are NaN in the output.

- **units**: units attribute if the data have none (required by xsdba methods).

- **return_model**: also return the corrector fitted on all years.

- **options**: futher options for xsbda methods (e.g. pooling, window_size, nquantiles, kind, interp).  See XsdbaCorrector for further information on all xsdba methods.

