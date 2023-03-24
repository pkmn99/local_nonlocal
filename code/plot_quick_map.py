import xarray as xr
from scipy import stats
import matplotlib.pyplot as plt

import cartopy.io.shapereader as shpreader
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import ShapelyFeature

#exp0='F2000climo_ctl'
#exp='F2000climo_Allgrass'
#exp='F2000climo_Allgrass_checkerboard'
#exp='F2000climo_Allgrass_minimal'
#var_group='clm2' #'clm2','cam'

def my_ttest(df1,df2):
    return stats.ttest_ind(df1,df2, equal_var = False)

# Vectorized t-test
def my_ttest_3d(df1,df2,dim='time'):
    f,p=xr.apply_ufunc(my_ttest, df1,df2, input_core_dims=[[dim],[dim]],
               output_core_dims=[[],[]],vectorize=True)
    return f,p

def load_data(exp,exp0='F2000climo_ctl',var_group='clm2',time_scale='yearmean'):
    # Load annual model output
    d0y=xr.open_dataset('../data/outputdata/%s.%s.h0.0001-0035-%s.nc'%(exp0,var_group,time_scale))
    d1y=xr.open_dataset('../data/outputdata/%s.%s.h0.0001-0035-%s.nc'%(exp,var_group,time_scale))
    
    # # # Seasonal 
    # # d0s=xr.open_dataset('../output_2.1.3/35years/F2000climo_f19_g16.%s.h0.0001-0035-seasmean.nc'%var_group)
    # # d1s=xr.open_dataset('../output_2.1.3/35years/%s.%s.h0.0001-0035-seasmean.nc'%(exp,var_group))

    # Calculate total precipitation  
    if var_group=='cam':
        d1y['PREC']=d1y['PRECC']+d1y['PRECL']
        d0y['PREC']=d0y['PRECC']+d0y['PRECL']

    print('Experiment %s data loaded'%exp)
    return d0y,d1y


# Input data is dataarray variable da0 from control and d1 from exp with time dimension 
    # Use shaply interface https://stackoverflow.com/questions/20990381/how-to-add-custom-shapefile-to-map-using-cartopy
    # with statistical values
def plot_change_map(da0,da1, ax,minmax=[],sig=False):
    # Load geographical data
    # china_shp=shpreader.Reader('../../data/China_gis/administration/China_southsea.shp')
    # china_shp=shpreader.Reader('../../data/China_gis/administration/CN-sheng-A-test.shp')
    world_shp=shpreader.Reader('../../data/China_gis/ne_110m_admin_0_countries/ne_110m_admin_0_countries.shp')

    world_feature = ShapelyFeature(world_shp.geometries(),
                                    ccrs.PlateCarree(), facecolor='none',edgecolor='grey', linewidth=0.5)
    # Land mask
    dl = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')

    ax.add_feature(world_feature)


#    ax.set_extent([60, 140, 10, 55], ccrs.Geodetic())
    
    #(da1-da0).mean('time').where(dl.LANDFRAC_PFT.values>0.5).plot(ax=ax)
    if minmax:
        (da1-da0).mean('time').where(dl.LANDFRAC_PFT.values>0.5).plot.contourf(ax=ax, levels=20,
                                                                               vmin=minmax[0],
                                                                               vmax=minmax[1],
                                                                              cmap='RdBu_r')
    else:
        (da1-da0).mean('time').where(dl.LANDFRAC_PFT.values>0.5).plot.contourf(ax=ax, levels=20,
                                                                              cmap='RdBu_r')
    # If conduct significant test, plot the significant grid box at 0.90 level with hatched line
    if sig:
        p=my_ttest_3d(da1, da0)[1]
        # Use hatch to represent significance
        cs = ax.contourf(p.lon,p.lat, 
                         (da1-da0).mean('time').where(dl.LANDFRAC_PFT.values>0.5).where(p<0.1),
                          colors='none',
                          hatches=['///'],
                          extend='lower',zorder=15
                          )

def make_plot(exp,var,time_scale='year',var_group='clm2',exp0='F2000climo_ctl'):
    if time_scale in ['DJF','MAM','JJA','SON']:
        d0y,d1y=load_data(exp,time_scale='seas',var_group=var_group,exp0=exp0)
    else:
        d0y,d1y=load_data(exp,time_scale=time_scale,var_group=var_group,exp0=exp0)

    # Plot change
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_axes([0, 0, 1, 1], projection=ccrs.PlateCarree(),
                      frameon=False)
    if time_scale in ['DJF','MAM','JJA','SON']:
        plot_change_map(d0y[var][d0y.time.dt.season==time_scale], 
                        d1y[var][d1y.time.dt.season==time_scale],
                        ax, minmax=[])
    else:
        plot_change_map(d0y[var], d1y[var], ax, minmax=[])
    ax.set_title(var)
    
    ax.set_global()
    plt.savefig('../figure/fig_%s_%s_%s_35years.png'%(exp,var,time_scale),dpi=300)
    print('figure saved')

if __name__ == "__main__":
#    make_plot('F2000climo_Allgrass','PREC',time_scale='JJA',var_group='cam')
     make_plot('F2000climo_Allgrass','TSA',var_group='clm2') # Global deforestation impact
     make_plot('F2000climo_Allgrass_checkerboard','TSA',var_group='clm2') # Global deforestation impact
