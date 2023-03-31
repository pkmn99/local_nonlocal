import xarray as xr
import numpy as np
"""
Calcualte local and nonlocal effect based on subgrid simulation
"""

# Calculate weigted mean based on pft fractions
def weighted_mean(da, w):
    d = da.copy().fillna(0)
    dw = (d.values * w.values).sum(axis=0)/w.values.sum(axis=0)
    dw[dw==0]=np.nan
    return  xr.DataArray(dw, coords=[da.lat, da.lon],
                         dims=["lat","lon"],name=da.name, attrs=da.attrs)

def save_data(var='TSA'):
    time_scale='monthly'
    d0 = xr.open_dataset('../data/outputdata/F2000climo_ctl.clm2.h1.0001-0035-%s-2d.nc'%time_scale)
    d1 = xr.open_dataset('../data/outputdata/F2000climo_Allgrass_minimal_tree.clm2.h1.0001-0035-%s-2d.nc'%time_scale)
    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    
    if var=='TSc':
        B = 5.670367*10**-8 # Stefan-Boltzmann constant
        d0y=((d0['FIRE']/(1*B))**0.25).mean(dim='time')
        d1y=((d1['FIRE']/(1*B))**0.25).mean(dim='time')
    else:
        d0y=d0[var].mean(dim='time')
        d1y=d1[var].mean(dim='time')

    # Local impact from subgrid differences: pft 1:11 forest+shrub, pft12:15 grass
    f =  weighted_mean(d1y[1:12], lc.PCT_NAT_PFT[1:12])#forest indef run 
    nf =  weighted_mean(d1y[12:15], lc.PCT_NAT_PFT[12:15])#nonforest in def run 
    f0 =  weighted_mean(d0y[1:12], lc.PCT_NAT_PFT[1:12]) #forest in control run
    
    d_l=nf-f # local impact
    d_nl=f-f0 # nonlocal impact
    
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])
    if var=='TSc':
        d_final[var+'_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                       'units':'K'}
        d_final[var+'_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                          'units':'K'}
    else:
        d_final[var+'_local'].attrs = d0[var].attrs
        d_final[var+'_nonlocal'].attrs = d0[var].attrs
    d_final.to_netcdf('../data/result/subgrid.%s.local_nonlocal.nc'%var)
    print('data saved')

if __name__=="__main__":
    save_data(var='TSc')
