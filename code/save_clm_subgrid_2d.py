import xarray as xr
import numpy as np

# Convert and save the 1d subgrid output of CLM to 2d data array
def convert_to_2d(df, var):
    a = np.full([df.time.shape[0],17, df.lat.shape[0], df.lon.shape[0]], np.nan)
    for t in range(df.time.shape[0]):
        for i in range(1,17):
            temp = np.full_like(df.landmask, np.nan)
            ind = df.pfts1d_itype_veg==i
            temp[df.pfts1d_jxy[ind]-1, df.pfts1d_ixy[ind]-1] = df[var][t,ind]
            a[t,i,:,:]=temp
    return  xr.DataArray(a, 
                     coords=[df.time, np.int32(list(range(0,17))), df.lat, df.lon],
                     dims=["time", "pft", "lat","lon"],name=var, attrs=df[var].attrs)

# Combine and save 2d array of different variables to a single dataframe 
def combine_dataframe(df, var_list=['EFLX_LH_TOT', 'FGR', 'FIRE', 'FSA', 'FSH', 'FSR', 'TLAI', 'TSA']):
    da_all=[convert_to_2d(df, i) for i in var_list]
    return  xr.merge(da_all)

def save_to_2d(fn_in, fn_out):
    dc2=xr.open_dataset(fn_in)
    dc2.load()
    df_final = combine_dataframe(dc2)
    df_final.to_netcdf(fn_out)
    print('file saved')

if __name__=="__main__":
#    exp='F2000climo_ctl'
    exp='F2000climo_Allgrass_minimal_tree'
    time_scale='monthly'
    fn_in='../data/outputdata/%s.clm2.h1.0001-0035-%s.nc'%(exp,time_scale)
    fn_out='../data/outputdata/%s.clm2.h1.0001-0035-%s-2d.nc'%(exp,time_scale)
    save_to_2d(fn_in, fn_out)
