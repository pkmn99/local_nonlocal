import xarray as xr
import numpy as np
"""
Calcualte local and nonlocal effect based on subgrid simulation
"""

# Calculate weigted mean across based on pft fractions
def weighted_mean(da, w):
    d = da.copy().fillna(0)
    dw = (d.values * w.values).sum(axis=0)/w.values.sum(axis=0)
    dw[dw==0]=np.nan
    return  xr.DataArray(dw, coords=[da.lat, da.lon],
                         dims=["lat","lon"],name=da.name, attrs=da.attrs)

def save_data(var='TSA'):
    d0 = xr.open_dataset('../data/outputdata/F2000climo_ctl.clm2.h1.0001-0035-yearmean-2d.nc')
    d1 = xr.open_dataset('../data/outputdata/F2000climo_Allgrass_minimal_tree.clm2.h1.0001-0035-yearmean-2d.nc')
    lc = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    
    d0y=d0[var].mean(dim='time')
    d1y=d1[var].mean(dim='time')

    # Local impact from subgrid differences: pft 1:11 forest+shrub, pft12:15 grass
    f =  weighted_mean(d1y[1:11], lc.PCT_NAT_PFT[1:11])
    nf =  weighted_mean(d1y[11:15], lc.PCT_NAT_PFT[11:15])
    f0 =  weighted_mean(d0y[1:11], lc.PCT_NAT_PFT[1:11])
    
    d_l=nf-f # local impact
    d_nl=f-f0 # nonlocal impact
    
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])
    d_final[var+'_local'].attrs = d0[var].attrs
    d_final[var+'_nonlocal'].attrs = d0[var].attrs
    d_final.to_netcdf('../data/result/subgrid.%s.local_nonlocal.nc'%var)
    print('data saved')

if __name__=="__main__":
    save_data(var='TSA')

