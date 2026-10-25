# SPDX-FileCopyrightText: 2024 European Centre for Medium-Range Weather Forecasts (ECMWF)
# SPDX-License-Identifier: Apache-2.0

# script that can downloads forecasted MJO indices from S2S.aux (dependent on Vitart)
import xarray as xr
import numpy as np
import pandas as pd
import ftplib
import os
import huracanpy as hpy
from acacia_s2s_toolkit import argument_output, download_S2Stc_tracks
from datetime import datetime, timedelta
from pathlib import Path
import yaml
import tarfile
import requests

# this is for sphinx - only functions listed here will have entries in readthedocs API
#__all__ = ["download_forecast_TCtracks", "download_reforecast_TCtracks"]

def open_aux_session():
    session = ftplib.FTP('aux.ecmwf.int')
    user, pwrd = download_S2Stc_tracks.get_s2sFTP_tokens()
    session.login(user=user,passwd=pwrd)
    return session

def download_fc_MJO_txt(model,fcdate):
    session = open_aux_session()
    session.cwd("RMMS")

    # get origin_id for filepath name. Need model for directory
    origin_id = argument_output.output_originID(model, fcdate)

    # create filepath
    mn = model.lower().strip() # model name
    yy = fcdate[:4] # year component
    mm = fcdate[4:6] # month component
    MJO_fn = f"z_s2s_rmm_{origin_id}_prod_rt_{fcdate}00.txt"

    filepath = f"{mn}/real-time/{yy}/{mm}/{MJO_fn}" # creates full filename

    local_filename = MJO_fn

    # retrieve MJO file
    with open(local_filename,'wb') as f:
        session.retrbinary(f"RETR {filepath}", f.write)

    print(f"File '{filepath}' has been downloaded.")
    session.quit()
    return local_filename # return the string used to save the file

def download_rfc_MJO_txt(model,fcdate,rfdate):
    session = open_aux_session()
    session.cwd("RMMS")

    # get origin_id for filepath name. Need model for directory
    origin_id = argument_output.output_originID(model, fcdate)

    # create filepath
    mn = model.lower().strip() # model name
    yy = fcdate[:4] # year component
    mm = fcdate[4:6] # month component
    yyrf = rfdate[:4]
    mmrf = rfdate[4:6]
    MJO_fn = f"z_s2s_rmm_{origin_id}_prod_rf_{rfdate}00_{fcdate}.txt"

    filepath = f"{mn}/reforecasts/{yy}/{mm}/{MJO_fn}" # creates full filename

    local_filename = MJO_fn

    # retrieve MJO file
    with open(local_filename,'wb') as f:
        session.retrbinary(f"RETR {filepath}", f.write)

    print(f"File '{filepath}' has been downloaded.")
    session.quit()
    return local_filename # return the string used to save the file

def single_MJO_fc_download(model,fcdate,leadtimes):
    # get origin_id
    origin_id = argument_output.output_originID(model, fcdate)

    # download single MJO forecast txt file
    MJO_fn = download_fc_MJO_txt(model,fcdate)

    # read downloaded MJO txt file.
    MJO_data = pd.read_csv(MJO_fn,skiprows=1,names=['lt','exptyp','member','RMM1','RMM2','Amplitude','Phase'],sep='\s+')
    # remove ensemble mean from text file (LATER OPTION, ENABLE DOWNLOAD OF SOLELY ENSEMBLE MEAN?). 
    MJO_data = MJO_data[MJO_data.exptyp != 'em']
    # set index based on lead time and ensemble member. convert to xarray.
    MJO_data = MJO_data.set_index(['lt','member']).to_xarray()
    # for ECMWF, 100th member appears as '**' (in text file as well).
    if origin_id == 'ecmf':
        MJO_data = MJO_data.assign_coords(member=['100' if m == '**' else m for m in MJO_data.member.values])
    # use integers for ensemble member rather than string, and sort by ensemble number. Remove member index of 0 by adding one to values
    MJO_data = (MJO_data.assign_coords(member_num=('member', MJO_data.member.astype(int).data+1)).swap_dims({'member': 'member_num'}).drop_vars('member').rename({'member_num': 'member'}).sortby('member'))
    # drop exptyp in coordinate
    MJO_data = MJO_data.drop_vars(['exptyp'])
    # compute time coordinate and replace lead time
    MJO_data = MJO_data.assign_coords(time=('lt',[datetime.strptime(fcdate,'%Y%m%d')+timedelta(hours=float(time)) for time in MJO_data.lt.data])).swap_dims({'lt':'time'}).drop_vars('lt').sortby('time')
    return MJO_data

def single_MJO_rfc_download(model,fcdate,rfdate,leadtimes):
    # get origin_id
    origin_id = argument_output.output_originID(model, fcdate)

    # download single MJO forecast txt file
    MJO_fn = download_rfc_MJO_txt(model,fcdate,rfdate)

    # read downloaded MJO txt file.
    MJO_data = pd.read_csv(MJO_fn,skiprows=1,names=['lt','exptyp','member','RMM1','RMM2','Amplitude','Phase'],sep='\s+')
    # remove ensemble mean from text file (LATER OPTION, ENABLE DOWNLOAD OF SOLELY ENSEMBLE MEAN?).
    MJO_data = MJO_data[MJO_data.exptyp != 'em']
    # set index based on lead time and ensemble member. convert to xarray.
    MJO_data = MJO_data.set_index(['lt','member']).to_xarray()
    # for ECMWF, 100th member appears as '**' (in text file as well).
    if origin_id == 'ecmf':
        MJO_data = MJO_data.assign_coords(member=['100' if m == '**' else m for m in MJO_data.member.values])
    # use integers for ensemble member rather than string, and sort by ensemble number. Remove member index of 0 by adding one to values
    MJO_data = (MJO_data.assign_coords(member_num=('member', MJO_data.member.astype(int).data+1)).swap_dims({'member': 'member_num'}).drop_vars('member').rename({'member_num': 'member'}).sortby('member'))
    # drop exptyp in coordinate
    MJO_data = MJO_data.drop_vars(['exptyp'])
    # compute time coordinate and replace lead time
    MJO_data = MJO_data.assign_coords(time=('lt',[datetime.strptime(fcdate,'%Y%m%d')+timedelta(hours=float(time)) for time in MJO_data.lt.data])).swap_dims({'lt':'time'}).drop_vars('lt').sortby('time')
    return MJO_data

def download_forecast_MJO(fcdate,model,origin_id,leadtime_hour,filename_save,fc_enslags):
    lag_i = 0
    all_fcs = []
    member_start = 1
    for lag in np.atleast_1d(fc_enslags):
        lag = int(lag)
        print (lag)
        # get correct leadtimes given selection of forecasts + lag
        leadtimes, convert_fcdate = argument_output.output_formatted_leadtimes(leadtime_hour,fcdate,'MJO',origin_id,lag=lag,fc_enslags=fc_enslags)
    
        print (leadtimes)
    
        date_obj = datetime.strptime(convert_fcdate, "%Y-%m-%d")
        convert_fcdate = date_obj.strftime("%Y%m%d")
    
        # download MJO forecast given lagged fc date and filtered leadtimes
        MJO_fc = single_MJO_fc_download(model,convert_fcdate,leadtimes) # job to do, figure out what to do with local_destination field
    
        nmem = MJO_fc.sizes["member"]
        MJO_fc = MJO_fc.assign_coords(member=np.arange(member_start,member_start+nmem))
        all_fcs.append(MJO_fc)
        member_start += nmem
    
    combined_allens = xr.concat(all_fcs,dim='member',join='outer')
    # detect whether any nan values are present, if so, remove those times
    bad_times = combined_allens.time[combined_allens.to_array().isnull().any(("variable", "member"))]
    if bad_times.size > 0:
        print (f'the following times {bad_times.values} contain NaN values. Most likely due to lagged ensemble use. Removing these times from forecast.')
        valid_times=~combined_allens.to_array().isnull().any(("variable", "member"))
        combined_allens = combined_allens.sel(time=valid_times)
    combined_allens.to_netcdf(f"{filename_save}.nc")
    print (f'Saved MJO indices in {filename_save}')
    # need to delete downloaded txt files!

def download_reforecast_MJO(fcdate,model,origin_id,leadtime_hour,filename_save,rf_enslags,rf_years,fc_time):
    lag_i = 0
    all_fcs = []
    member_start = 1
    for lag in np.atleast_1d(rf_enslags):
        lag = int(lag)
        print (lag)
        # get correct leadtimes given selection of forecasts + lag
        # convert fcdate
        lagged_fcdate = datetime.strptime(fcdate, '%Y%m%d')+timedelta(days=lag) # work out what the lagged fcdate is.
        convert_fcdate = lagged_fcdate.strftime('%Y%m%d') # convert that date to YYYYMMDD format

        rf_model_date, rfyears = argument_output.check_and_output_all_hc_arguments('MJO',origin_id,convert_fcdate,rf_years) # get the reforecast model date version plus the number of reforecast years

        all_rfyrs = []

        for rfyear in np.atleast_1d(sorted(rfyears)):
            new_rfdate = f"{rfyear}{convert_fcdate[4:8]}"
            
            # download MJO forecast given lagged fc date and filtered leadtimes
            MJO_rfc = single_MJO_rfc_download(model,convert_fcdate,new_rfdate,leadtimes) # job to do, figure out what to do with local_destination field

            nmem = MJO_rfc.sizes["member"]
            MJO_rfc = MJO_rfc.assign_coords(member=np.arange(member_start,member_start+nmem))
            all_rfyrs.append(MJO_rfc)
        combined_allens = xr.concat(all_rfyrs,dim='year',join='outer')
        all_fcs.append(combined_allens)
        member_start += nmem

    combined_allens = xr.concat(all_fcs,dim='member',join='outer')
    # detect whether any nan values are present, if so, remove those times
    bad_times = combined_allens.time[combined_allens.to_array().isnull().any(("variable", "member"))]
    if bad_times.size > 0:
        print (f'the following times {bad_times.values} contain NaN values. Most likely due to lagged ensemble use. Removing these times from forecast.')
        valid_times=~combined_allens.to_array().isnull().any(("variable", "member"))
        combined_allens = combined_allens.sel(time=valid_times)
    combined_allens.to_netcdf(f"{filename_save}.nc")
    print (f'Saved reforecast MJO indices in {filename_save}')
    # need to delete downloaded txt files!

