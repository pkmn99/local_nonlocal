import xarray as xr
from plot_quick_map import load_data
from local_subgrid import get_total_effect

# calculate local and nonlocal impact from offline simulation
def save_data(var='TSc'):
    [c, g] = load_data(exp='I2000Clm50Sp_Allgrass',exp0='I2000Clm50Sp_ctl',var_group='clm2',time_scale='monthly')
    
    B = 5.670367*10**-8 # Stefan-Boltzmann constant
    c['TSc']=((c['FIRE']/(1*B))**0.25)
    g['TSc']=((g['FIRE']/(1*B))**0.25)
    
    d_l = (g.TSc-c.TSc)
    d_total = get_total_effect()
    d_nl=d_total-d_l # nonlocal impact as residue
    d_final = xr.merge([d_l.rename(var+'_local'),d_nl.rename(var+'_nonlocal')])

    if var=='TSc':
        d_final[var+'_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                       'units':'K'}
        d_final[var+'_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                          'units':'K'}
    else:
        d_final[var+'_local'].attrs = d0[var].attrs
        d_final[var+'_nonlocal'].attrs = d0[var].attrs

    d_final.to_netcdf('../data/result/offline.%s.local_nonlocal.mon.nc'%var)
    print('data saved')

if __name__=="__main__":
    save_data(var='TSc')
