import xarray as xr
import numpy as np

# Convert and save the 1d subgrid output of CLM to 2d data array
def convert_to_2d(df, var):
    a = np.full([df.time.shape[0],17, df.lat.shape[0], df.lon.shape[0]], np.nan)
    for t in range(df.time.shape[0]):
        for i in range(0,17):
            temp = np.full_like(df.landmask, np.nan)
            ind = df.pfts1d_itype_veg==i
            temp[df.pfts1d_jxy[ind]-1, df.pfts1d_ixy[ind]-1] = df[var][t,ind]
            a[t,i,:,:]=temp
    return  xr.DataArray(a, 
                     coords=[df.time, np.int32(list(range(0,17))), df.lat, df.lon],
                     dims=["time", "pft", "lat","lon"],name=var, attrs=df[var].attrs)

# Modify for 350x720 hourly data with seperate mapping and data file 04/21/2021
def convert_to_2d_360(df, var, df_au):
    a = np.full([df.time.shape[0],17, df_au.lat.shape[0], df_au.lon.shape[0]], np.nan)
    for t in range(df.time.shape[0]):
        for i in range(0,17):
            temp = np.full_like(df_au.landmask, np.nan)
            ind = df_au.pfts1d_itype_veg==i
            temp[df_au.pfts1d_jxy[ind]-1, df_au.pfts1d_ixy[ind]-1] = df[var][t,ind]
            a[t,i,:,:]=temp
    return  xr.DataArray(a, 
                     coords=[df.time, np.int32(list(range(0,17))), df_au.lat, df_au.lon],
                     dims=["time", "pft", "lat","lon"],name=var, attrs=df[var].attrs)


# Combine and save 2d array of different variables to a single dataframe 
def combine_dataframe(df, var_list):
    da_all=[convert_to_2d(df, i) for i in var_list]
    return  xr.merge(da_all)

# MOdify for 360x720 houraly data
def combine_dataframe_360(df, var_list):
    au=xr.open_dataset('../data/outputdata/auxi.nc')
    da_all=[convert_to_2d_360(df, i, au) for i in var_list]
    return  xr.merge(da_all)

def save_to_2d(fn_in, fn_out,var_list=['EFLX_LH_TOT', 'FGR', 'FIRE', 'FSA', 'FSH', 'FSR', 'TLAI', 'TSA']):
    dc2=xr.open_dataset(fn_in)
    dc2.load()
    df_final = combine_dataframe(dc2,var_list=var_list)
    df_final.to_netcdf(fn_out)
    print('file saved')

def save_to_2d_360(fn_in, fn_out,var_list=['EFLX_LH_TOT','FSR', 'FSH']):
    dc2=xr.open_dataset(fn_in)
    dc2.load()
    df_final = combine_dataframe_360(dc2,var_list=var_list)
    df_final.to_netcdf(fn_out)
    print('file saved %s'%fn_out)

if __name__=="__main__":
    exp='I2000Clm50Sp_360x720cru'
#    I2000Clm50Sp_360x720cru.clm2.h1.0015-12-31-00000.nc
    time_scale='monthly'
#    fn_in='../data/outputdata/%s.clm2.h1.0001-0035-%s.nc'%(exp,time_scale)
#    fn_out='../data/outputdata/%s.clm2.h1.0001-0035-%s-2d.nc'%(exp,time_scale)
#    fn_in='../data/outputdata/bowen_data/I2000Clm50Sp_360x720cru.clm2.h1.0006.monmean_day.nc'
#    fn_in='../data/outputdata/bowen_data/test_FSH_mean.nc'
#    fn_out='../data/outputdata/bowen_data/test_FSH_mean_2d.nc'
    fn_in='../data/outputdata/bowen_data/bowen_ratio_25_75_mean_1d.nc'
    fn_out='../data/outputdata/bowen_data/bowen_ratio_25_75_mean_2d.nc'
    save_to_2d_360(fn_in, fn_out,var_list=['FSH','EFLX_LH_TOT','BOWEN_RATIO'])
    fn_in='../data/outputdata/bowen_data/bowen_ratio_25_75_median_1d.nc'
    fn_out='../data/outputdata/bowen_data/bowen_ratio_25_75_median_2d.nc'
    save_to_2d_360(fn_in, fn_out,var_list=['FSH','EFLX_LH_TOT','BOWEN_RATIO'])
