import xarray as xr
import numpy as np
from plot_quick_map import load_data

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
    return d_local_interp, d_nonlocal_interp

def save_data(var='TSA',exp='F2000climo_Allgrass_checkerboard',var_group='clm2'):
    # Load necessary data
    d0y, d1y = load_data(exp,exp0='F2000climo_ctl',var_group=var_group)
    d = (d1y[var] - d0y[var]).mean(dim='time')
    
    # Get the deforestation grid box
    lc0 = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr-Allgrass_checkerboard.nc')
    def_loc=(lc0-lc).PCT_NAT_PFT[1:12,:,:].sum(dim='natpft')>0 # True means deforestation grids

    d_l,d_nl=local_nonlocal_interp(d, def_loc)
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])
    d_final[var+'_local'].attrs = d0y[var].attrs
    d_final[var+'_nonlocal'].attrs = d0y[var].attrs
#    d_final.to_netcdf('../data/result/%s.%s.local_nonlocal.nc'%(exp,var))
    d_final.to_netcdf('../data/result/interpolation.%s.local_nonlocal.nc'%var)
    print('data saved')

if __name__=='__main__':
    save_data(var='CLDTOT',var_group='cam')
