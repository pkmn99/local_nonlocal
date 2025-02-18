import xarray as xr
import numpy as np
from plot_quick_map import load_data
from local_nonlocal_ibm import load_def_frac

# Use interpolation method to separate local and nonlocal effects from checkerboard pattern
def local_nonlocal_interp(d,def_loc):
    landmask = xr.open_dataset('../data/outputdata/I2000Clm50Sp_Allgrass_minimal_tree.clm2.h1.0001-0035-yearmean.nc')['landmask']
    d_nonlocal=d.copy()
    d_nonlocal=xr.where(def_loc.values,np.nan,d_nonlocal)
    # Use mean of both directions (bilinear? in my case)
    d_nonlocal_interp=((d_nonlocal.interpolate_na(dim="lat", method="linear")
                        +d_nonlocal.interpolate_na(dim="lon", method="linear"))/2)
#     d_nonlocal_interp =xr.where(np.isnan(d),np.nan,d_nonlocal_interp) # Mask ocean values during intepolation
    d_nonlocal_interp =xr.where(landmask.values!=1,np.nan,d_nonlocal_interp) # Mask ocean values during intepolation

    d_local=(d-d_nonlocal_interp)
    d_local=xr.where(~def_loc.values,np.nan,d_local)
    d_local_interp=((d_local.interpolate_na(dim="lat", method="linear")
                     +d_local.interpolate_na(dim="lon", method="linear"))/2)
#     d_local_interp =xr.where(np.isnan(d),np.nan,d_local_interp)
    d_local_interp =xr.where(landmask.values!=1,np.nan,d_local_interp)
    return d_local_interp.transpose("time",...), d_nonlocal_interp.transpose("time",...)

def save_data(var='TSA',exp0='F2000climo_ctl',exp1='F2000climo_Allgrass_checkerboard',var_group='clm2',time_scale='monthly'):
    # Load necessary data
    d0y, d1y = load_data(exp1,exp0=exp0,var_group=var_group,time_scale=time_scale)
    if var=='TSc':
        B = 5.670367*10**-8 # Stefan-Boltzmann constant
        d=((d1y['FIRE']/(1*B))**0.25-(d0y['FIRE']/(1*B))**0.25)
    else:
        d = (d1y[var] - d0y[var])
    
    # Get the deforestation grid box
#    lc0 = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
#    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr-Allgrass_checkerboard.nc')
#    def_loc=(lc0-lc).PCT_NAT_PFT[1:12,:,:].sum(dim='natpft')>0 # True means deforestation grids
    def_loc = load_def_frac(exp0=exp0,exp1=exp1)
    d_l,d_nl=local_nonlocal_interp(d, def_loc>0)
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])
    if var=='TSc':
        d_final[var+'_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                       'units':'K'}
        d_final[var+'_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                          'units':'K'}
    else:
        d_final[var+'_local'].attrs = d0y[var].attrs
        d_final[var+'_nonlocal'].attrs = d0y[var].attrs
#    d_final.to_netcdf('../data/result/%s.%s.local_nonlocal.nc'%(exp,var))
    d_final.sel(time=slice('0006','0035')).to_netcdf('../data/result/interpolation.%s.local_nonlocal.mon.nc'%var)
    print('data saved')

if __name__=='__main__':
#    save_data(var='CLDTOT',var_group='cam')
    save_data(var='TSc',var_group='clm2')
