import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from local_nonlocal_ibm import load_data
from local_subgrid import get_total_effect
# Use Tair of the bottom layer of atmosphtere as nonlocal effect
def save_data(exp0='F2000climo_ctl',exp1='F2000climo_Allgrass'):

    # Total actual surface temperature difference
    d_total = get_total_effect()

    g_cam = load_data(exp1,var_group='cam',time_scale='monthly')
    c_cam = load_data(exp0,var_group='cam',time_scale='monthly')
    
    g_clm = load_data(exp1,var_group='clm2',time_scale='monthly')
    c_clm = load_data(exp0,var_group='clm2',time_scale='monthly')
    
    # Use atmoshperic air temperature change as nonlocal impact Ta
    tair_c=c_cam['T'][:,31] # air temperature of the bottom atmosphere level
    tair_g=g_cam['T'][:,31] # air temperature of the bottom atmosphere level
    
    tair_c=tair_c.drop_vars('lev')
    tair_g=tair_g.drop_vars('lev')
    tair_c['lat']=d_total.lat
    tair_g['lat']=d_total.lat
    
    nonlocal_tair=tair_g-tair_c # changes in backgroud cliamte by Ta 
    
    local=d_total-nonlocal_tair

    d_final = xr.merge([local.rename('TSc_local'),nonlocal_tair.rename('TSc_nonlocal')])
    
    d_final['TSc_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                  'units':'K'}
    d_final['TSc_nonlocal'].attrs = {'long name': 'derived from changes in bottom atmos tair',
                                     'units':'K'}
    d_final.sel(time=slice('0006','0035')).to_netcdf('../data/result/Tair.TSc.local_nonlocal.mon.nc')
    print('data saved')

if __name__=="__main__":
    save_data()
