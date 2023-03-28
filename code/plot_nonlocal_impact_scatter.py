import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from plot_local_impact_scatter import my_scatter

#plot scatterplot comparision of nonlocal impact 
def make_plot():
    # load data for comparison
    di=xr.open_dataset('../data/result/interpolation.TSc.local_nonlocal.nc')# checkerboard
    ds=xr.open_dataset('../data/result/subgrid.TSc.local_nonlocal.nc')#copuled subgrid
    dibm=xr.open_dataset('../data/result/ibm.TSc.local_nonlocal.nc')# IBM
    doff=xr.open_dataset('../data/result/offline.TSc.local.nc')# offline
    dsoff=xr.open_dataset('../data/result/subgrid.TSc.local.nc')# subgrid offline

    var='TSc_nonlocal'
    # order is 1 ds_r, 2 dsoff_r, 3 di, 4 doff, 5 dibm, 6 ibm-Tair

    fig = plt.figure(figsize=(14,11))
    ax11 = fig.add_subplot(5, 5, 1)
    my_scatter(dsoff[var].values.flatten(), ds[var].values.flatten(),ax11,x0=0.2)
    
    ax21 = fig.add_subplot(5, 5, 6)
    my_scatter(di[var].values.flatten(), ds[var].values.flatten(),ax21,x0=0.2)
    
    ax31 = fig.add_subplot(5, 5, 11)
    my_scatter(doff[var].values.flatten(), ds[var].values.flatten(),ax31,x0=0.2)
    
    ax41 = fig.add_subplot(5, 5, 16)
    my_scatter(dibm[var].values.flatten(), ds[var].values.flatten(),ax41,x0=0.2)
    
    ax51 = fig.add_subplot(5, 5, 21)
    my_scatter(dibm['Tair_nonlocal'].values.flatten(), ds[var].values.flatten(),ax51,x0=0.2)
    
    ax22 = fig.add_subplot(5, 5, 7)
    my_scatter(di[var].values.flatten(), dsoff[var].values.flatten(),ax22,x0=0.2)
    
    ax32 = fig.add_subplot(5, 5, 12)
    my_scatter(doff[var].values.flatten(), dsoff[var].values.flatten(),ax32,x0=0.2)
    
    ax42 = fig.add_subplot(5, 5, 17)
    my_scatter(dibm[var].values.flatten(), dsoff[var].values.flatten(),ax42,x0=0.2)
    
    ax52 = fig.add_subplot(5, 5, 22)
    my_scatter(dibm['Tair_nonlocal'].values.flatten(), dsoff[var].values.flatten(),ax52,x0=0.2)
    
    ax33 = fig.add_subplot(5, 5, 13)
    my_scatter(doff[var].values.flatten(), di[var].values.flatten(),ax33,x0=0.2)
    
    ax43 = fig.add_subplot(5, 5, 18)
    my_scatter(dibm[var].values.flatten(), di[var].values.flatten(),ax43,x0=0.2)
    
    ax53 = fig.add_subplot(5, 5, 23)
    my_scatter(dibm['Tair_nonlocal'].values.flatten(), di[var].values.flatten(),ax53,x0=0.2)
    
    ax44 = fig.add_subplot(5, 5, 19)
    my_scatter(dibm[var].values.flatten(), doff[var].values.flatten(),ax44,x0=0.2)
    
    ax54 = fig.add_subplot(5, 5, 24)
    my_scatter(dibm['Tair_nonlocal'].values.flatten(), doff[var].values.flatten(),ax54,x0=0.2)
    
    ax55 = fig.add_subplot(5, 5, 25)
    my_scatter(dibm['Tair_nonlocal'].values.flatten(), dibm[var].values.flatten(),ax55,x0=0.2)
    
    # set range etc
    for a in [ax11, ax21, ax31,ax41, ax21, ax22, ax31,ax32,ax33,
              ax41,ax42,ax43,ax44,ax51,ax52,ax53,ax54,ax55]:
        a.set_yticks(np.arange(-4,3.1,2))
        a.set_xticks(np.arange(-4,3.1,2))
        a.set_xlim(-4, 3)
        a.set_ylim(-4, 3)
        a.plot([-4,3],[-4,3],'--b',lw=1)
    
    # add label  
    [a.set_ylabel('Coupled subgrid ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax11,ax21,ax31,ax41,ax51]]
    [a.set_ylabel('Offline subgrid ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax22,ax32,ax42,ax52]]
    [a.set_ylabel('Checkerboard ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax33,ax43,ax53]]
    [a.set_ylabel('Offline ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax44,ax54]]
    ax55.set_ylabel('IBM-Tair ($^\circ$C)', fontsize=10,labelpad=0)
    
    [a.set_xlabel('IBM-Tair ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax51,ax52,ax53,ax54,ax55]]
    [a.set_xlabel('IBM ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax41,ax42,ax43,ax44]]
    [a.set_xlabel('Offline ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax31,ax32,ax33]]
    [a.set_xlabel('Checkerboard ($^\circ$C)', fontsize=10,labelpad=0) for a in [ax21,ax22]]
    ax11.set_xlabel('Offline subgrid ($^\circ$C)', fontsize=10,labelpad=0)
    
    fig.subplots_adjust(hspace=0.35, wspace=0.4)
    
    plt.savefig('../figure/nonlocal_impact_scatter20230328.png',dpi=300,bbox_inches='tight')
    print('fig saved')

if __name__=="__main__":
    make_plot()
