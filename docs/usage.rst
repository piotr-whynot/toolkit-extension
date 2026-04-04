Usage
=====

Basic Example
-------------

Here's a simple example of using toolkit-extension::

    import toolkit_extension as te

    #define work directory

    data_dir="./data"

    nominal_date="2026-03-01"
    target_domain="madagascar"
    fcst_var="tp"
    fcst_model="ECMWF"
    obs_file="./data/PRCPTOT_day_CHC_CHIRPS-2.0-0p25_merged.nc"
    obs_var="PRCPTOT"  #name of the variable stored in the obs_file
    domain_shape_file="./data/madagascar.geojson" # a file with vector data showing the boundaries of the domain - will be used in plotting only


    agg_window="5D"  #data will be aggregated into blocks of this size
    agg_method="sum" #and it will be a sum of daily values over these blocks

    # this defines regridding parameters for matching the obs and forecast grids. 
    # it assumes obs grid is finer, so fine_to_coarse will give forecast grid, coarse_to_fine will give obs grid
    # raise_if_missing will stop processing if the overlap between obs and forecast grids of less than fractional 
    # threshold
    grid_alignment_kwargs=dict(
        direction="fine_to_coarse",
        method="conservative",
        raise_if_missing=True,
        threshold=0.9
    )

    # this defines time alignment parameter for matching the obs and hindcast time series. 
    # raise_if_missing will stop processing if any of the hindcast days are not covered by observations 

    time_alignment_kwargs={"raise_if_missing":True}


    #calling processing function
    hindcast_lt,forecast_lt,obs_lt=te.preprocess_forecast(
        nominal_date=nominal_date,
        download_dir=data_dir,
        target_domain=target_domain,
        fcst_var=fcst_var,
        fcst_model=fcst_model,
        agg_window=agg_window,
        agg_method=agg_method,
        obs_file=obs_file,
        obs_var=obs_var,
        gridalignment_kwargs=grid_alignment_kwargs,
        timealignment_kwargs=time_alignment_kwargs,
        verbose=True
    )


More Details
------------

Some more explanation here...
