# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# Check that requested variables are appropriate and are compatiable with WMO lead centre.
from acacia_s2s_toolkit.variable_dict import s2s_variables
from difflib import get_close_matches
from datetime import datetime, timedelta
from acacia_s2s_toolkit import argument_output
import numpy as np

def check_requested_variable(variable):
    '''check requested variable matches abbreviations used for the S2S database.
    No return - ECDS variable outputted in variable_output.py
    '''
    
    # Flatten all variables from nested dictionary
    all_vars = []
    for category_dict in s2s_variables.values():
        for subcategory_vars in category_dict.values():
            all_vars.extend(subcategory_vars)

    if variable in all_vars:
        return True  # Variable is valid
    else:
        # Try suggesting closest matches
        suggestions = get_close_matches(variable, all_vars, n=3)
        suggestion_msg = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
        
        raise ValueError(
            f"Invalid variable: '{variable}' is not in the S2S database.{suggestion_msg}"
        )

def check_model_name(model,fcdate):
    # Flatten all models from nested dictionary
    df = argument_output.read_lookup_table(fcdate)
    all_models = list(df['Model'].values)

    if model in all_models:
        return True  # Variable is valid
    else:
        # Try suggesting closest matches
        suggestions = get_close_matches(model, all_models, n=3)
        suggestion_msg = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""

        raise ValueError(
            f"Invalid model: '{model}' is not currently compatible with this toolbox. {suggestion_msg}"
            )

def check_fcdate(fcdate,origin_id):
    """
    Check if fcdate is a valid date string in the format 'YYYYMMDD'.
    Returns a datetime.date object if valid, or raises an error if invalid.
    """
    if not isinstance(fcdate, str):
        raise ValueError(f"[ERROR] Forecast date must be a string, got {type(fcdate)}.")
        return None

    try:
        date_obj = datetime.strptime(fcdate, '%Y%m%d')
    except ValueError:
        raise ValueError(f"[ERROR] '{fcdate}' is not in the correct format 'YYYYMMDD'.")
        return None

    # given originID check requested fcdate is permitted
    now = datetime.utcnow()

    # get time difference between requested time and now. 
    time_diff = now - date_obj # how many days bigger is now compared to requested time.

    # using origin name and origin_latency_hours dictionary, check time_diff is larger than numbers of hours.
    # get latency period
    min_numhours = argument_output.get_single_parameter(origin_id,fcdate,'Delay')

    # check min_numhours is smaller than time_diff. for instant, you cannot request an ECMWF forecast after 24 hours. 
    if time_diff < timedelta(hours=float(min_numhours)):
       raise ValueError(f"[ERROR] The time difference between now and requested forecast date {time_diff} is smaller than the required minimum amount of time for originID {origin_id}.") 

    # check that requested forecast date, matches avaliable forecast initilisations
    # get weekday forecast initialisations
    weekdays_aval = argument_output.get_single_parameter(origin_id,fcdate,'fcFreq')

    # check that the forecast date weekday is avaliable for that model
    fcdate_weekday = date_obj.weekday()+1 # Monday = 1 etc..
    if fcdate_weekday not in weekdays_aval:
        raise ValueError(f"[ERROR] The chosen forecast initialisation date is not avaliable for the chosen model. Origin ID: {origin_id}")

def check_dataformat(data_format):
    if data_format not in ['grib','netcdf']:
        raise ValueError(f"[ERROR] The chosen data format is not avaliable. Please use 'grib' or 'netcdf'")

def check_leadtime_hours(leadtime_hour,variable,origin_id,fcdate):
    # is the maximum lead time smaller or equal to forecat end time
    end_time = argument_output.get_single_parameter(origin_id,fcdate,'fcLength')

    if np.max(leadtime_hour) > end_time:
        raise ValueError(f"[ERROR] You are requesting a leadtime greater than end of forecast, {end_time} hours")

    if np.min(leadtime_hour) < 0:
        raise ValueError(f"[ERROR] You are requesting a negative leadtime")

    # the check depends on time resolution
    time_resolution = argument_output.get_timeresolution(variable)

    if time_resolution.endswith('6hrly'):
        output_freq = 6
    else:
        output_freq = 24

    if not np.all(leadtime_hour % output_freq == 0):
        raise ValueError(f"[ERROR] You are requesting a leadtime that is not compatible with output frequency. Output frequency of the desired variable is {output_freq} hours. You are requesting the following {leadtime_hour}.")

def check_plevs(plevs,variable):
    # first get maximum plevs avaliable
    max_plevs = argument_output.output_plevs(variable)

    # Check all requested plevs are valid
    if not np.all(np.isin(plevs, max_plevs)):
        invalid = [p for p in plevs if p not in max_plevs]
        raise ValueError(
            f"[ERROR] Invalid pressure level(s) requested: {invalid}. "
            f"Available levels for '{variable}' are: {max_plevs}."
        )

def check_area_selection(area):
    # Go through each component of the area (N, W, S, E)
    if area[0] < -90 or area[0] > 90:
        raise ValueError(
            f"[ERROR] Invalid northern latitude '{area[0]}'. Must be between -90 and 90."
            )

    if area[1] < -180 or area[1] > 180:
        raise ValueError(
            f"[ERROR] Invalid western longitude '{area[1]}'. Must be between -180 and 180."
            )

    if area[2] < -90 or area[2] > 90:
        raise ValueError(
            f"[ERROR] Invalid southern latitude '{area[2]}'. Must be between -90 and 90."
            )

    if area[3] < -180 or area[3] > 180:
        raise ValueError(
           f"[ERROR] Invalid eastern longitude '{area[3]}'. Must be between -180 and 180."
        )

    # is north > south?
    if area[0] < area[2]:
        raise ValueError(
                f"[ERROR] Northern latitude {area[0]} must be greater than southern latitude {area[2]}."
            )
    if area[1] > area[3]:
         raise ValueError(
                f"[ERROR] Western longitude {area[1]} must be smaller than eastern longitude {area[3]}."
            )

def check_fc_enslags(fc_enslags):
    """
    Check that all values in fc_enslags are non-positive integers (i.e., ≤ 0 and whole numbers).
    Raises ValueError if check fails.
    """

    # Case 1: a single integer
    if isinstance(fc_enslags, int):
        if fc_enslags > 0:
            raise ValueError("All ensemble lags (fc_enslags) must be integers ≤ 0 (e.g., [0, -1, -2]).")

    # Case 2: iterable of integers, a list
    else:
        if not all(isinstance(lag, int) and lag <= 0 for lag in fc_enslags):
            raise ValueError("All ensemble lags (fc_enslags) must be integers ≤ 0 (e.g., [0, -1, -2]).")

def check_requested_reforecast_years(rf_years,origin_id,fc_date):
    ''' 
    Check that the requested reforecast years are able to download
    '''
    # first get full set of reforecast years
    full_rf_years = argument_output.get_hindcast_year_span(origin_id,fc_date)
   
    # check all years in rf_years are in full_rf_years
    if not all(year in full_rf_years for year in rf_years):
        raise ValueError(f"All requested reforecast years {rf_years} are not avaliable. Avaliable years are {full_rf_years}.")
