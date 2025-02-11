import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
# from scipy import stats
# from matplotlib.colors import ListedColormap, Normalize
from matplotlib.colors import Normalize
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def make_plot():
    # load data for comparison
   # di=xr.open_dataset('../data/result/interpolation.TSc.local_nonlocal.nc')# checkerboard
    di=xr.open_dataset('../data/result/subgrid.TSc.local_nonlocal.nc')#copuled subgrid
   # ds=xr.open_dataset('../data/result/subgrid.TSc.local_nonlocal.nc')#copuled subgrid
    ds=xr.open_dataset('../data/result/subgrid_r.TSc.local_nonlocal.nc')#copuled subgrid
    dibm=xr.open_dataset('../data/result/ibm.TSc.local_nonlocal.nc')# IBM
    dibm_r=xr.open_dataset('../data/result/ibm_r.TSc.local_nonlocal.nc')# IBM with reversed ctl and lc exp
    doff=xr.open_dataset('../data/result/offline.TSc.local.nc')# offline
    dsoff=xr.open_dataset('../data/result/subgrid.TSc.local.nc')# subgrid offline
    
    # Get the deforestation grid box
    lc0 = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
    
    # forest change fraction
    def_loc=lc0.PCT_NAT_PFT[1:12,:,:].sum(dim='natpft')
    
    # multiple local impact by fractions for subgrid results
#    ds_r=ds * def_loc.values / 100
    ds_r=ds# * def_loc.values / 100
    dsoff_r=dsoff * def_loc.values / 100
    
    var='TSc_local'

    # begin plot 
    levels=np.arange(-2,2.1,0.2)
    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(3, 2, 1, projection=ccrs.PlateCarree())
    ds_r[var].plot(ax=ax1,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)
    
    ax2 = fig.add_subplot(3, 2, 2, projection=ccrs.PlateCarree())
    dsoff_r[var].plot(ax=ax2,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)
    
    ax3 = fig.add_subplot(3, 2, 3, projection=ccrs.PlateCarree())
    di[var].plot(ax=ax3,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)
    
    ax4 = fig.add_subplot(3, 2, 4, projection=ccrs.PlateCarree())
    doff[var].plot(ax=ax4,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)
    
    ax5 = fig.add_subplot(3, 2, 5, projection=ccrs.PlateCarree())
    dibm[var].plot(ax=ax5,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)

    ax6 = fig.add_subplot(3, 2, 6, projection=ccrs.PlateCarree())
    (-dibm_r[var]).plot(ax=ax6,cmap='RdBu_r',vmin=-2,vmax=2,add_colorbar=False,levels=levels)
    
    # set boundary and display extent 
    [a.coastlines() for a in [ax1,ax2,ax3,ax4,ax5,ax6]]
    [a.set_extent([-180, 180, -60, 80]) for a in [ax1,ax2,ax3,ax4,ax5,ax6]]
    
    ax1.set_title('Coupled subgrid')
    ax2.set_title('Offline subgrid')
    ax3.set_title('Checkerboard')
    ax4.set_title('Offline')
    ax5.set_title('IBM')
    ax6.set_title('IBM_reverse')
    
    # Add colorbar
    cbar1_pos = [ax5.get_position().x0, ax5.get_position().y0-0.04,
                 ax5.get_position().width*1.0, 0.015]
    cax1 = fig.add_axes(cbar1_pos)
    cb1 = mpl.colorbar.ColorbarBase(ax=cax1, cmap='RdBu_r', norm=Normalize(vmin=-2, vmax=2),
                                    orientation='horizontal', ticks=np.arange(-2, 2.1, 0.5)) #cmap=plt.get_cmap('hot')
    cb1.set_label('$\Delta$TS ($^\circ$C)', fontsize=12)
    
    # Add panel label
    ax1.text(-0.06, 1.01, 'a', fontsize=14, transform=ax1.transAxes, fontweight='bold')
    ax2.text(-0.06, 1.01, 'b', fontsize=14, transform=ax2.transAxes, fontweight='bold')
    ax3.text(-0.06, 1.01, 'c', fontsize=14, transform=ax3.transAxes, fontweight='bold')
    ax4.text(-0.06, 1.01, 'd', fontsize=14, transform=ax4.transAxes, fontweight='bold')
    ax5.text(-0.06, 1.01, 'e', fontsize=14, transform=ax5.transAxes, fontweight='bold')
    
    plt.savefig('../figure/local_impact_comparision20250126.png',dpi=300,bbox_inches='tight')
    print('fig saved')

if __name__=="__main__":
    make_plot()
