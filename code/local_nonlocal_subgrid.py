import xarray as xr
import numpy as np
from local_nonlocal_ibm import load_def_frac
"""
Calcualte local and nonlocal effect based on subgrid simulation
"""

# Calculate weigted mean based on pft fraction*s
def weighted_mean_2d(da, w):
    d = da.copy().fillna(0)
    dw = (d.values * w.values).sum(axis=0)/w.values.sum(axis=0)
    dw[dw==0]=np.nan
    return  xr.DataArray(dw, coords=[da.lat, da.lon],
                         dims=["lat","lon"],name=da.name, attrs=da.attrs)

# Calculate weigted mean based on pft fraction*s
# for 3d dim with time dimension
def weighted_mean_3d(da, w):
    d = da.copy().fillna(0)
    dw = (d.values * w.values).sum(axis=1)/w.values.sum(axis=0)
    dw[dw==0]=np.nan
    return  xr.DataArray(dw, coords=[da.time, da.lat, da.lon],
                         dims=["time","lat","lon"],name=da.name, attrs=da.attrs)

def save_data(exp0='F2000climo_ctl',exp1='F2000climo_Allgrass_minimal_tree',var='TSA'):
    time_scale='monthly'
    d0 = xr.open_dataset('../data/outputdata/%s.clm2.h1.0001-0035-%s-2d.nc'%(exp0,time_scale))
    d1 = xr.open_dataset('../data/outputdata/%s.clm2.h1.0001-0035-%s-2d.nc'%(exp1,time_scale))
    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    # forest change fraction
#    def_frac=lc.PCT_NAT_PFT[1:12,:,:].sum(dim='natpft')
    def_frac=load_def_frac(exp0=exp0,exp1=exp1)
    
    if var=='TSc':
        B = 5.670367*10**-8 # Stefan-Boltzmann constant
        d0y=((d0['FIRE']/(1*B))**0.25)#.mean(dim='time')
        d1y=((d1['FIRE']/(1*B))**0.25)#.mean(dim='time')
    else:
        d0y=d0[var]#.mean(dim='time')
        d1y=d1[var]#.mean(dim='time')

    # Local impact from subgrid differences: pft 1:11 forest+shrub, pft12:15 grass
    f =  weighted_mean_3d(d1y[:,1:12], lc.PCT_NAT_PFT[1:12])#forest in def run 
    nf =  weighted_mean_3d(d1y[:,12:15], lc.PCT_NAT_PFT[12:15])#nonforest in def run 
    f0 =  weighted_mean_3d(d0y[:,1:12], lc.PCT_NAT_PFT[1:12]) #forest in control run
    
    d_l=(nf-f) # local impact
    d_lw=(nf-f)*def_frac.values/100 # local impact, mutiply by deforestation fraction
    d_nl=f-f0 # nonlocal impact
    
    d_final = xr.merge([d_lw.rename(var+'_local'),d_l.rename(var+'_localuw'),d_nl.rename(var+'_nonlocal')])
    if var=='TSc':
        d_final[var+'_localuw'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                       'units':'K'}
        d_final[var+'_local'].attrs = {'long name': 'calculated surface temperature from FIRE multiply by deforestation fraction',
                                       'units':'K'}
        d_final[var+'_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                          'units':'K'}
    else:
        d_final[var+'_local'].attrs = d0[var].attrs
        d_final[var+'_nonlocal'].attrs = d0[var].attrs
    if exp0=='F2000climo_ctl': 
       d_final.sel(time=slice('0006','0035')).to_netcdf('../data/result/subgrid.%s.local_nonlocal.mon.nc'%var)
    if exp0=='F2000climo_Allgrass_minimal_tree':
       d_final.to_netcdf('../data/result/subgrid_r.%s.local_nonlocal.nc'%var)

    print('data saved')

if __name__=="__main__":
    save_data(exp0='F2000climo_ctl',exp1='F2000climo_Allgrass_minimal_tree',var='TSc')
