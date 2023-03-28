import xarray as xr
import numpy as np
from local_nonlocal_subgrid import weighted_mean
from local_nonlocal_ibm import load_data

"""
Calcualte local and nonlocal effect based on subgrid simulation
"""
# get total effect from the coupled simulation
def get_total_effect(var='TSc'):
    g_clm = load_data('F2000climo_Allgrass',var_group='clm2',time_scale='monthly')
    c_clm = load_data('F2000climo_ctl',var_group='clm2',time_scale='monthly')

    # Calculate Ts from outgoing longwave radiation
    B = 5.670367*10**-8 # Stefan-Boltzmann constant

    lwout_c=c_clm['FIRE'] # emitted longwave radiation
    lwout_g=g_clm['FIRE'] # emitted longwave radiation

    ts_c=(lwout_c/(1*B))**0.25
    ts_g=(lwout_g/(1*B))**0.25

    ts_diff = ts_g-ts_c # actual surface temperature difference
    return ts_diff

def save_data(var='TSc'):
    time_scale='monthly'
    d1 = xr.open_dataset('../data/outputdata/I2000Clm50Sp_Allgrass_minimal_tree.clm2.h1.0001-0035-%s-2d.nc'%time_scale)
    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    
    if var=='TSc':
        B = 5.670367*10**-8 # Stefan-Boltzmann constant
        d1y=((d1['FIRE']/(1*B))**0.25).mean(dim='time')
    else:
        d1y=d1[var].mean(dim='time')

    # Local impact from subgrid differences: pft 1:11 forest+shrub, pft12:15 grass
    f =  weighted_mean(d1y[1:12], lc.PCT_NAT_PFT[1:12])
    nf = weighted_mean(d1y[12:15], lc.PCT_NAT_PFT[12:15])
    
    d_l=nf-f # local impact
    d_total = get_total_effect()
    # nonlocal impact as residue, mutiple local impact with fraction before calculating nonlocal
    d_nl=d_total.mean(dim='time')-d_l*lc.PCT_NAT_PFT[1:12,:,:].sum(dim='natpft').values/100
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])

    if var=='TSc':
        d_final[var+'_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                       'units':'K'}
        d_final[var+'_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                          'units':'K'}
    else:
        d_final[var+'_local'].attrs = d0[var].attrs
        d_final[var+'_nonlocal'].attrs = d0[var].attrs
    d_final.to_netcdf('../data/result/subgrid.%s.local20230324.nc'%var)
    print('data saved')

if __name__=="__main__":
    save_data(var='TSc')
