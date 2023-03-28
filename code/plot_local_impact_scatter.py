import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Return values excluding nan values
def xy_no_nan(x,y):
    return x[~np.isnan(x)&~np.isnan(y)],y[~np.isnan(x)&~np.isnan(y)]

# Calculate correlation excluding nan values
def my_corr(x,y):
    x1,y1=xy_no_nan(x, y)
    r, p = pearsonr(x1, y1)
    return r,p
# calculate rmse excluding nan values
def my_rmse(x,y):
    x1,y1=xy_no_nan(x, y)
    return (((x1 - y1)**2).sum()/y1.shape[0])**0.5

# define scatterplot style
def my_scatter(x,y,ax, x0=0.3):
    ax.scatter(x, y, s=5,c='grey')
    r,p=my_corr(x, y)
    rmse_temp = my_rmse(x,y)
    ax.text(x0,0.025,'r=%.2f rmse=%.2f'%(r,rmse_temp), transform=ax.transAxes)

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

    var='TSc_local'
    # order is 1 ds_r, 2 dsoff_r, 3 di, 4 doff, 5 dibm

    fig = plt.figure(figsize=(13,10))
    ax11 = fig.add_subplot(4, 4, 1)
    my_scatter(dsoff_r[var].values.flatten(), ds_r[var].values.flatten(),ax11)
    
    ax21 = fig.add_subplot(4, 4, 5)
    my_scatter(di[var].values.flatten(), ds_r[var].values.flatten(),ax21)
    
    ax31 = fig.add_subplot(4, 4, 9)
    my_scatter(doff[var].values.flatten(), ds_r[var].values.flatten(),ax31)
    
    ax41 = fig.add_subplot(4, 4, 13)
    my_scatter(dibm[var].values.flatten(), ds_r[var].values.flatten(),ax41)
    
    ax22 = fig.add_subplot(4, 4, 6)
    my_scatter(di[var].values.flatten(), dsoff_r[var].values.flatten(),ax22)
    
    ax32 = fig.add_subplot(4, 4, 10)
    my_scatter(doff[var].values.flatten(), dsoff_r[var].values.flatten(),ax32)
    
    ax42 = fig.add_subplot(4, 4, 14)
    my_scatter(dibm[var].values.flatten(), dsoff_r[var].values.flatten(),ax42)
    
    ax33 = fig.add_subplot(4, 4, 11)
    my_scatter(doff[var].values.flatten(), di[var].values.flatten(),ax33)
    
    ax43 = fig.add_subplot(4, 4, 15)
    my_scatter(dibm[var].values.flatten(), di[var].values.flatten(),ax43)
    
    ax44 = fig.add_subplot(4, 4, 16)
    my_scatter(dibm[var].values.flatten(), doff[var].values.flatten(),ax44)
    
    # set range etc
    for a in [ax11, ax21, ax31,ax41, ax21, ax22, ax31,ax32,ax33,ax41,ax42,ax43,ax44]:
        a.set_yticks(np.arange(-1,2.6,1))
        a.set_xticks(np.arange(-1,2.6,1))
        a.set_xlim(-2, 3)
        a.set_ylim(-2, 3)
        a.plot([-2,3],[-2,3],'--b',lw=1)

    # add label 
    [a.set_ylabel('Coupled subgrid ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax11,ax21,ax31,ax41]]
    [a.set_ylabel('Offline subgrid ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax22,ax32,ax42]]
    [a.set_ylabel('Checkerboard ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax33,ax43]]
    ax44.set_ylabel('Offline ($^\circ$C)', fontsize=9,labelpad=0)
    
    [a.set_xlabel('IBM ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax41,ax42,ax43,ax44]]
    [a.set_xlabel('Offline ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax31,ax32,ax33]]
    [a.set_xlabel('Checkerboard ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax21,ax22]]
    ax11.set_xlabel('Offline subgrid ($^\circ$C)', fontsize=10,labelpad=0)
    
    fig.subplots_adjust(hspace=0.35, wspace=0.4) 

    plt.savefig('../figure/local_impact_scatter20230328.png',dpi=300,bbox_inches='tight')
    print('fig saved')

if __name__=="__main__":
    make_plot()
