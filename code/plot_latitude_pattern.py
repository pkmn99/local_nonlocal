import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def make_plot():
    # load data for comparison
    di=xr.open_dataset('../data/result/interpolation.TSc.local_nonlocal.nc')# checkerboard
    ds=xr.open_dataset('../data/result/subgrid.TSc.local_nonlocal.nc')#copuled subgrid
    dibm=xr.open_dataset('../data/result/ibm.TSc.local_nonlocal.nc')# IBM
    doff=xr.open_dataset('../data/result/offline.TSc.local.nc')# offline
    dsoff=xr.open_dataset('../data/result/subgrid.TSc.local.nc')# subgrid offline

    # Get the deforestation grid box
    lc0 = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')

    # forest change fraction
    def_loc=lc0.PCT_NAT_PFT[1:12,:,:].sum(dim='natpft')

    # multiple local impact by fractions for subgrid results
    ds_r=ds * def_loc.values / 100
    dsoff_r=dsoff * def_loc.values / 100

    fig, [ax1,ax2]=plt.subplots(1,2,figsize=(12,4))

    var='TSc_local'
    ds_r[var].mean(dim='lon').plot(ax=ax1)
    dsoff_r[var].mean(dim='lon').plot(ax=ax1)
    di[var].mean(dim='lon').plot(ax=ax1)
    doff[var].mean(dim='lon').plot(ax=ax1)
    dibm[var].mean(dim='lon').plot(ax=ax1)
    ax1.plot([-60,80],[0,0],'k--',lw=0.5)
    ax1.set_xlim([-60,80])
    ax1.set_ylim([-1.7,1.5])
    ax1.set_title('Local impact')
    ax1.set_ylabel('Temperature change ($^\circ$C)')
    ax1.set_xlabel('Latitude')
    ax1.legend(['Coupled subgrid','Offline subgrid','Checkerboard','Offline','IBM'])

    var='TSc_nonlocal'
    ds[var].mean(dim='lon').plot(ax=ax2)
    dsoff[var].mean(dim='lon').plot(ax=ax2)
    di[var].mean(dim='lon').plot(ax=ax2)
    doff[var].mean(dim='lon').plot(ax=ax2)
    dibm[var].mean(dim='lon').plot(ax=ax2)
    dibm['Tair_nonlocal'].mean(dim='lon').plot(ax=ax2)
    ax2.plot([-60,80],[0,0],'k--',lw=0.5)
    ax2.set_xlim([-60,80])
    ax2.set_ylim([-1.7,1.5])
    ax2.set_title('Nonlocal impact')
    ax2.set_ylabel('Temperature change ($^\circ$C)')
    ax2.set_xlabel('Latitude')
    ax2.legend(['Coupled subgrid','Offline subgrid','Checkerboard','Offline','IBM','IBM_Tair'])

    ax1.text(-0.06, 1.05, 'a', fontsize=14, transform=ax1.transAxes, fontweight='bold')
    ax2.text(-0.06, 1.05, 'b', fontsize=14, transform=ax2.transAxes, fontweight='bold')
    
    plt.savefig('../figure/latitude_pattern20230331.png',dpi=300,bbox_inches='tight')
    print('fig saved')

if __name__=="__main__":
    make_plot()
