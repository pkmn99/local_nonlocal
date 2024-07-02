import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

def load_data(exp,var_group='clm2',time_scale='monthly'):
    d=xr.open_dataset('../data/outputdata/%s.%s.h0.0001-0035-%s.nc'%(exp,var_group,time_scale))
    return d

def save_data():
    g_cam = load_data('F2000climo_Allgrass',var_group='cam',time_scale='monthly')
    c_cam = load_data('F2000climo_ctl',var_group='cam',time_scale='monthly')
    
    g_clm = load_data('F2000climo_Allgrass',var_group='clm2',time_scale='monthly')
    c_clm = load_data('F2000climo_ctl',var_group='clm2',time_scale='monthly')
    
    # Calculate Ts from outgoing longwave radiation
    B = 5.670367*10**-8 # Stefan-Boltzmann constant
    
    h_c = c_clm['FSH'] # sensible heat 
    h_g = g_clm['FSH'] # sensible heat
    
    le_c = c_clm['EFLX_LH_TOT'] # latent heat
    le_g = g_clm['EFLX_LH_TOT'] # latent heat
    
    g_c=c_clm['FGR'] # ground flux
    g_g=g_clm['FGR'] # ground flux
    
    lwout_c=c_clm['FIRE'] # emitted longwave radiation
    lwout_g=g_clm['FIRE'] # emitted longwave radiation
    
    lwin_c=c_clm['FLDS'] # incoming longwave radiation
    lwin_g=g_clm['FLDS'] # incoming longwave radiation
    
    swin_c=c_clm['FSDS'] # incoming shortwave radiation
    swin_g=g_clm['FSDS'] # incoming shortwave radiation
    
    swout_c=c_clm['FSR'] # reflected shortwave radiation
    swout_g=g_clm['FSR'] # reflected shortwave radiation
    
    # calculate surface temperature from longwave radiation
    ts_c=(lwout_c/(1*B))**0.25
    ts_g=(lwout_g/(1*B))**0.25
    
    # calculate abledo
    alb_c=(c_clm.FSR/c_clm.FSDS)
    alb_g=(g_clm.FSR/g_clm.FSDS)

    # atmoshperic air temperature as Ta
    tair_c=c_cam['T'][:,31] # air temperature of the bottom atmosphere level
    tair_g=g_cam['T'][:,31] # air temperature of the bottom atmosphere level
    
    tair_c=tair_c.drop_vars('lev')
    tair_g=tair_g.drop_vars('lev')
    tair_c['lat']=ts_c.lat
    tair_g['lat']=ts_g.lat
    
    ts_diff = ts_g-ts_c # actual surface temperature difference
    
    # Calculate necessary variables for IBM method
    ra_c =(ts_c-tair_c)/h_c * 1.2250 * 1003.5 # aerodynamic roughness
    bo_c = h_c/le_c # bowen ratio
    lam_c = 1/(4*B*tair_c**3) #lambda term
    F_c = (1.2250 * 1003.5 *lam_c)/ra_c*(1+1/bo_c) # energy redistribution factor
    
    ra_g =(ts_g-tair_g)/h_g * 1.2250 * 1003.5
    bo_g = h_g/le_g
    lam_g = 1/(4*B*tair_g**3)
    F_g = (1.2250 * 1003.5 *lam_g)*(1+1/bo_g)/ra_g
    
    # Net radiation
    Rn_c = swin_c - swout_c + lwin_c - 1*B*tair_c**4
    Rn_g = swin_g - swout_g + lwin_g - 1*B*tair_g**4
    
    # Use CTL to calculate IBM
    ra = ra_c
    bo = bo_c
    F = F_c
    Rn=Rn_c
    
    delta_f1= -lam_c *1.2250*1003.5 *(1+1/bo) * (ra_g-ra_c)/ra**2
    delta_f2= -lam_c *1.2250*1003.5 * (bo_g-bo_c)/(ra*bo**2)
    
    # Decomposed tempreature
    dts1 = lam_c /(1+F) * (c_clm.FSDS*-(alb_g-alb_c))
    
    dts2 = -lam_c/((1+F)**2) * Rn * delta_f1
    
    dts3 = -lam_c/(1+F)**2 * Rn * delta_f2
    
    dts4=tair_g-tair_c # changes in backgroud cliamte by Ta 
    
    dts=dts1+dts2+dts3

    local=dts.where((dts>ts_diff.quantile(0.001))&(dts<ts_diff.quantile(0.999))).mean(dim='time') # remove outliers
    nonlocal_tair=dts4.mean(dim='time')
    nonlocal_tsc=ts_diff.mean(dim='time')-local
    
    d_final = xr.merge([local.rename('TSc_local'),nonlocal_tsc.rename('TSc_nonlocal'),nonlocal_tair.rename('Tair_nonlocal')])
    
    d_final['TSc_local'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                  'units':'K'}
    d_final['TSc_nonlocal'].attrs = {'long name': 'calculated surface temperature from FIRE',
                                     'units':'K'}
    d_final['Tair_nonlocal'].attrs = tair_c.attrs
    d_final.to_netcdf('../data/result/ibm.TSc.local_nonlocal.nc')
    print('data saved')

if __name__=="__main__":
    save_data()
