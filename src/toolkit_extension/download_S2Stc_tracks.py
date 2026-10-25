# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# script that can downloaded forecasted S2S TC tracks
import xarray as xr
import numpy as np
import pandas as pd
import ftplib
import os
import huracanpy as hpy
from acacia_s2s_toolkit import argument_output
from datetime import datetime, timedelta
from pathlib import Path
import yaml
import tarfile
import requests

# this is for sphinx - only functions listed here will have entries in readthedocs API
__all__ = ["download_forecast_TCtracks", "download_reforecast_TCtracks"]

def safe_extract_tar(archive_path, extract_dir="."):
    archive_path = Path(archive_path)
    extract_dir = Path(extract_dir).resolve()

    with tarfile.open(archive_path) as tar:
        for member in tar.getmembers():
            member_path = (extract_dir / member.name).resolve()
            if not str(member_path).startswith(str(extract_dir) + os.sep) and member_path != extract_dir:
                raise ValueError(f"Unsafe path in tar archive: {member.name}")
        tar.extractall(path=extract_dir)

def empty_member_dataset():
    return xr.Dataset(
        data_vars={
            'track_id': ('record', np.array([], dtype='int64')),
            'time': ('record', np.array([], dtype='datetime64[ns]')),
            'lat': ('record', np.array([], dtype='float64')),
            'lon': ('record', np.array([], dtype='float64')),
            'wind': ('record', np.array([], dtype='float64')),
            'pres': ('record', np.array([], dtype='float64')),
            'lat_wind': ('record', np.array([], dtype='float64')),
            'lon_wind': ('record', np.array([], dtype='float64')),
            'basin': ('record', np.array([], dtype='U3')),
        },
        coords={
            'ens_mem': ('record', np.array([], dtype='int64')),
        }
    )

def get_s2sFTP_tokens():
    # open text file contain aux site passwords
    url = 'https://sites.ecmwf.int/ecbox/acacia_s2s_toolkit/s/api/v2/pub/files/S2Saux_user_pwrd'
    r = requests.get(url)
    r.raise_for_status()
    creds = r.json()

    return creds["user"], creds["password"]

def load_fc_tracks(model,fcdate):
    # load in data
    session = ftplib.FTP('aux.ecmwf.int')
    user, pwrd = get_s2sFTP_tokens()
    session.login(user=user,passwd=pwrd)
    session.cwd("TCYC")
    
    # create filepath
    mn = model.lower().strip() # model name
    yy = fcdate[:4] # year component
    mm = fcdate[4:6] # month component
    TC_fn = f"TC.{fcdate}"
    
    filepath = f"{mn}/real-time/{yy}/{mm}/{TC_fn}" # creates full filename
    local_filename = TC_fn
    
    # retrieve the full year file
    with open(local_filename,'wb') as f:
        session.retrbinary(f"RETR {filepath}", f.write)
    
    print(f"File '{filepath}' has been downloaded.")
    
    session.quit()
    return local_filename

def load_rffc_tracks(model,rfdate,fcdate):
    # load in data
    session = ftplib.FTP('aux.ecmwf.int')
    user, pwrd = get_s2sFTP_tokens()
    session.login(user=user,passwd=pwrd)
    session.cwd("TCYC")

    # create filepath
    mn = model.lower().strip() # model name
    yy_fc = fcdate[:4] # year component
    mm_fc = fcdate[4:6] # month component
    yy_rf = rfdate[:4] # yr component 
    mm_rf = rfdate[4:6] # month component of rfdate
    
    TC_fn = f"TC.{model}.{rfdate}.{fcdate}" # filename based on rfdate and fcdate

    filepath = f"{mn}/reforecasts/{yy_fc}/{mm_fc}/{TC_fn}" # creates full filename
    local_filename = TC_fn

    # retrieve the full year file
    with open(local_filename,'wb') as f:
        session.retrbinary(f"RETR {filepath}", f.write)

    print(f"File '{filepath}' has been downloaded.")

    session.quit()
    return local_filename

def download_reforecast_TCtracks(fcdate,model,origin_id,leadtime_hour,filename_save,rf_enslags,rf_years,fc_time):
    '''
    Downloads tracked TCs within reforecasts
    '''
    # extract multiple ens lags
    lag_i = 0
    all_fcs = []
    for lag in np.atleast_1d(sorted(rf_enslags)):
        lag = int(lag)
        # convert fcdate
        lagged_fcdate = datetime.strptime(fcdate, '%Y%m%d')+timedelta(days=lag) # work out what the lagged fcdate is.
        convert_fcdate = lagged_fcdate.strftime('%Y%m%d') # convert that date to YYYYMMDD format

        rf_model_date, rfyears = argument_output.check_and_output_all_hc_arguments('TC_TRACKS',origin_id,convert_fcdate,rf_years) # get the reforecast model date version plus the number of reforecast years

        all_rfyrs = []

        # loop across rf_years
        for rfyear in np.atleast_1d(sorted(rfyears)):
            new_rfdate = f"{rfyear}{convert_fcdate[4:8]}" # label yyyymmdd
            fn = load_rffc_tracks(model,new_rfdate,rf_model_date) # load in forecast TC track file from Frederic's FTP site

            # untar the file. will give all ensemble members. for ECMWF 101.
            safe_extract_tar(fn)

            short_name = origin_id  # need to get short name. add to look-up file
            num_pert_hcs = argument_output.get_single_parameter(origin_id,convert_fcdate,'rfNumEns')+1

            # loop through all basins
            basins = ['atl','aus','cnp','enp','nin','sin','spc','wnp']

            for num in range(num_pert_hcs):
                # untar first ens mem
                safe_extract_tar(f"{short_name}.{convert_fcdate}.{new_rfdate}.{num}")

                all_storms = []

                # loop through basins and read files
                for basin in basins:
                    # open up single basin file
                    with open(basin, "r") as f:
                        lines = f.readlines()

                    if len(lines) > 1:
                        storms = hpy.load(basin,source='ecmwf')
                        if basin == 'sin' or basin == 'spc' or basin == 'aus':
                            storms['lat'] = storms['lat']*-1.0 # convert to SH latitudes
                        n = storms.sizes['record']
                        # add dimension with basin
                        storms = storms.assign_coords(basin=('record', np.repeat(basin, n)))

                        all_storms.append(storms)

                # concatenate along the existing 'record' axis
                if not all_storms: # np storms exist so skip this member
                    combined = xr.Dataset().assign_coords(record=[])
                else:
                    combined = xr.concat(
                        all_storms,
                        dim='record',
                        join='outer',                       # allow variables that may be missing in some basins
                        combine_attrs='drop_conflicts'      # avoid attribute conflicts
                    )

                if 'basin' in combined.coords:
                    combined = combined.reset_coords('basin')

                # enforce a consistent variable ordering
                desired_vars = ['track_id','time','lat','lon','wind','pres','lat_wind','lon_wind','basin']
                available = [v for v in desired_vars if v in combined.variables]
                combined = combined[available]

                # Skip if empty
                if 'record' not in combined.sizes or combined.sizes['record'] == 0:
                    # empty dataset for this ensemble member
                    combined = empty_member_dataset()
                    all_fcs.append(combined)
                    continue

                n = combined.sizes['record']
                lagged_ens_num = num+(num_pert_hcs*lag_i)+1 # create a new ensemble number based on lag*num_ens_mem
                combined = combined.assign_coords(ens_mem=('record', np.repeat(lagged_ens_num, n)))

                # select TCs within chosen leadtimes
                # using requested leadtime_hour and fcdate, only keep storms on requested valid times
                t0 = datetime.strptime(new_rfdate, "%Y%m%d")
                valid_dates = [t0 + timedelta(hours=int(h)) for h in leadtime_hour]

                # make same version
                valid_np = np.array(valid_dates, dtype='datetime64[ns]')
                times_np = combined['time'].values.astype('datetime64[ns]')

                # mask
                mask = np.isin(times_np, valid_np)

                # Subset  
                combined = combined.isel(record=mask)
                all_rfyrs.append(combined)

                # remove extracted file
                os.remove(f"{short_name}.{convert_fcdate}.{new_rfdate}.{num}")
            os.remove(f"{fn}") # plus remove original tar file downloaded
        # concentnate all rfs
        combined_allrfs = xr.concat(
                  all_rfyrs,
            dim='record',
            join='outer',                       # allow variables that may be missing in some basins
            combine_attrs='drop_conflicts'      # avoid attribute conflicts
        )
        if fc_time==True:
            combined_allrfs['time'] = combined_allrfs['time']-pd.Timedelta(days=lag)

        all_fcs.append(combined_allrfs) # append to all_fcs array 

        lag_i = lag_i + 1

    # concatenate along the existing 'record' axis
    combined_allens = xr.concat(
            all_fcs,
            dim='record',
            join='outer',                       # allow variables that may be missing in some basins
            combine_attrs='drop_conflicts'      # avoid attribute conflicts
        )

    combined_allens = combined_allens.reset_coords('ens_mem')
    combined_allens = combined_allens[['track_id', 'time', 'lat', 'lon', 'wind', 'pres', 'lat_wind', 'lon_wind', 'basin','ens_mem']]

    vars_to_drop = ['lat_wind','lon_wind']
    combined_allens = combined_allens.drop_vars([v for v in vars_to_drop if v in combined_allens.variables])

    # using the combined_allens xarray, select TCs within chosen leadtimes
    print (combined_allens)

    combined_allens.to_netcdf(f"{filename_save}.nc")

    print(f"downloaded reforecast TC tracks for {num_pert_hcs*lag_i} ensemble members for rf_years {rfyears}, model {model} at fc_date {fcdate}. Saved as netcdf output under filename, {filename_save}")

    # then will need to remove basin files from last ensemble member
    for basin in basins:
        os.remove(f"{basin}")

def download_forecast_TCtracks(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags): 
    # extract multiple ens lags
    lag_i = 0
    all_fcs = []
    for lag in np.atleast_1d(sorted(fc_enslags)):
        lag = int(lag)
        print (lag)
        leadtimes, convert_fcdate = argument_output.output_formatted_leadtimes(leadtime_hour,fcdate,'TC_TRACKS',origin_id,lag=lag,fc_enslags=fc_enslags)

        print (leadtimes)

        date_obj = datetime.strptime(convert_fcdate, "%Y-%m-%d")
        convert_fcdate = date_obj.strftime("%Y%m%d")

        fn = load_fc_tracks(model,convert_fcdate) # load in forecast TC track file from Frederic's FTP site
        
        # untar the file. will give all ensemble members. for ECMWF 101.
        safe_extract_tar(fn)
        
        short_name = origin_id  # need to get short name. add to look-up file
        num_ens_mems = argument_output.get_single_parameter(origin_id,convert_fcdate,'fcNumEns')+1
        
        # loop through all basins
        basins = ['atl','aus','cnp','enp','nin','sin','spc','wnp']
        
        for num in range(num_ens_mems):
            # untar first ens mem
            safe_extract_tar(f"{short_name}.{convert_fcdate}.{num}")
            
            all_storms = []
            
            # loop through basins and read files
            for basin in basins:
                # open up single basin file
                with open(basin, "r") as f:
                    lines = f.readlines()
            
                if len(lines) > 1:
                    storms = hpy.load(basin,source='ecmwf')
                    if basin == 'sin' or basin == 'spc' or basin == 'aus':
                        storms['lat'] = storms['lat']*-1.0 # convert to SH latitudes
                    n = storms.sizes['record']
                    # add dimension with basin
                    storms = storms.assign_coords(basin=('record', np.repeat(basin, n)))
            
                    all_storms.append(storms)
            
            # concatenate along the existing 'record' axis
            if not all_storms: # np storms exist so skip this member
                combined = xr.Dataset().assign_coords(record=[])
            else:
                combined = xr.concat(
                    all_storms,
                    dim='record',
                    join='outer',                       # allow variables that may be missing in some basins
                    combine_attrs='drop_conflicts'      # avoid attribute conflicts
                )
            
            if 'basin' in combined.coords:
                combined = combined.reset_coords('basin')
            
            # enforce a consistent variable ordering
            desired_vars = ['track_id','time','lat','lon','wind','pres','lat_wind','lon_wind','basin']
            available = [v for v in desired_vars if v in combined.variables]
            combined = combined[available]
            
            # Skip if empty
            if 'record' not in combined.sizes or combined.sizes['record'] == 0:
                # empty dataset for this ensemble member
                combined = empty_member_dataset()
                all_fcs.append(combined)
                continue

            n = combined.sizes['record']
            lagged_ens_num = num+(num_ens_mems*lag_i)+1 # create a new ensemble number based on lag*num_ens_mem
            combined = combined.assign_coords(ens_mem=('record', np.repeat(lagged_ens_num, n)))
        
            all_fcs.append(combined)
        
            # remove extracted file
            os.remove(f"{short_name}.{convert_fcdate}.{num}")
        lag_i = lag_i + 1

    # concatenate along the existing 'record' axis
    combined_allens = xr.concat(
            all_fcs,
            dim='record',
            join='outer',                       # allow variables that may be missing in some basins
            combine_attrs='drop_conflicts'      # avoid attribute conflicts
        )
    
    combined_allens = combined_allens.reset_coords('ens_mem')
    combined_allens = combined_allens[['track_id', 'time', 'lat', 'lon', 'wind', 'pres', 'lat_wind', 'lon_wind', 'basin','ens_mem']]
    
    vars_to_drop = ['lat_wind','lon_wind']
    combined_allens = combined_allens.drop_vars([v for v in vars_to_drop if v in combined_allens.variables])

    # using the combined_allens xarray, select TCs within chosen leadtimes
    print (combined_allens)

    # using requested leadtime_hour and fcdate, only keep storms on requested valid times
    t0 = datetime.strptime(fcdate, "%Y%m%d")
    valid_dates = [t0 + timedelta(hours=int(h)) for h in leadtime_hour]

    # make same version
    valid_np = np.array(valid_dates, dtype='datetime64[ns]')
    times_np = combined_allens['time'].values.astype('datetime64[ns]')

    # mask
    mask = np.isin(times_np, valid_np)

    # Subset  
    combined_allens = combined_allens.isel(record=mask)

    combined_allens.to_netcdf(f"{filename_save}.nc")
    
    print(f"downloaded TC tracks for {num_ens_mems*lag_i} ensemble members for model {model} at fc_date {fcdate}. Saved as netcdf output under filename, {filename_save}")

    # then will need to remove basin files from last ensemble member
    for basin in basins:
        os.remove(f"{basin}")
    os.remove(f"{fn}") # plus remove original tar file downloaded
