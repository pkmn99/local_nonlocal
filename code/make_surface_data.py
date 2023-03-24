import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

"""
Create land surface variable for local and nonlocal CESM/CLM experiment
Modify the "PCT_NAT_PFT" variable to convert tree/shrub to grass pft 
"""

# Name could be "allgrass", "allgrass_checkerboard", and "allgrass_minimal_tree"
def create_land_cover(name):
    #d = xr.open_dataset('../../cesm/inputdata/inputdata_min/lnd/clm2/surfdata_map/release-clm5.0.18/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304.nc')
    d = xr.open_dataset('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr.nc')
     
    # create Allgrass world
    # sum of grass pft fraction
    grass = d.PCT_NAT_PFT[12:15].sum(dim='natpft')
    
    # nongrass pft fraction
    if name=='allgrass_minimal_tree':
        #  create minimal trees/shrub pft fractions based on Allgrass world
        # temp variable ngpft
        ngpft=d.PCT_NAT_PFT[1:12].copy()
        # change each tree/shrub pft to a minimal
        # if > 0.01 then 0.01; if <0.01, original fraction; this ensures that no negative pft values 
        for i in range(0,11):
            ngpft[i]=xr.where(ngpft[i]>0.01, 0.01, ngpft[i])
        # sum of nongrass pft fraction (trees and shrubs)
        nongrass = d.PCT_NAT_PFT[1:12].sum(dim='natpft')-ngpft.sum(dim='natpft')
    else:
        # sum of nongrass pft fraction (trees and shrubs)
        nongrass = d.PCT_NAT_PFT[1:12].sum(dim='natpft')

    # Assign nongrass pft fraction to grass based on their existing ratios
    # Get the modified grass pft fraction 
    pft12=d.PCT_NAT_PFT[12,:,:]*(nongrass/grass+1)
    pft13=d.PCT_NAT_PFT[13,:,:]*(nongrass/grass+1)
    pft14=d.PCT_NAT_PFT[14,:,:]*(nongrass/grass+1)

    if name=='allgrass': 
        # Change nongrass pft to zero
        d.PCT_NAT_PFT[1:12,:,:] = 0 
        
        # Assign mofided grass pft fraction to data
        d.PCT_NAT_PFT[12]=pft12
        d.PCT_NAT_PFT[13]=pft13
        d.PCT_NAT_PFT[14]=pft14
        
        # Fill na with zero and save the modified file
        d['PCT_NAT_PFT'] = d.PCT_NAT_PFT.fillna(0)
        d.to_netcdf('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr-Allgrass.nc',format='NETCDF3_64BIT')
        print("max pft is %f, and min pft is %f"%(d.PCT_NAT_PFT.max(),d.PCT_NAT_PFT.min()))
        print('land cover data for %s saved'%name)

    if name=='allgrass_checkerboard': 
        # Create checkboard moasic, using the method below
        # https://stackoverflow.com/questions/2169478/how-to-make-a-checkerboard-in-numpy
        checker_mask=np.indices([96,144]).sum(axis=0) % 2 
        
        # Change nongrass pft to zero based on the checkerboard pattern
        d.PCT_NAT_PFT[1:12,:,:] = xr.where(checker_mask==1, 0, d.PCT_NAT_PFT[1:12,:,:])
        
        # Assign modified grass pft fraction
        d.PCT_NAT_PFT[12,:,:] = xr.where(checker_mask==1, pft12, d.PCT_NAT_PFT[12,:,:])
        d.PCT_NAT_PFT[13,:,:] = xr.where(checker_mask==1, pft13, d.PCT_NAT_PFT[13,:,:])
        d.PCT_NAT_PFT[14,:,:] = xr.where(checker_mask==1, pft14, d.PCT_NAT_PFT[14,:,:])
        
        d['PCT_NAT_PFT'] = d.PCT_NAT_PFT.fillna(0)
        d.to_netcdf('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr-Allgrass_checkerboard.nc',format='NETCDF3_64BIT')
        print("max pft is %f, and min pft is %f"%(d.PCT_NAT_PFT.max(),d.PCT_NAT_PFT.min()))
        print('land cover data for %s saved'%name)
    
    if name=='allgrass_minimal_tree': 
        # Assign mofided grass pft fraction to data
        d.PCT_NAT_PFT[1:12]=ngpft
        d.PCT_NAT_PFT[12]=pft12
        d.PCT_NAT_PFT[13]=pft13
        d.PCT_NAT_PFT[14]=pft14
        d['PCT_NAT_PFT'] = d.PCT_NAT_PFT.fillna(0)
        d.to_netcdf('../data/inputdata/surfdata_1.9x2.5_hist_16pfts_Irrig_CMIP6_simyr2000_c190304_noirr-Allgrass_minimal_tree.nc',format='NETCDF3_64BIT')
        print("max pft is %f, and min pft is %f"%(d.PCT_NAT_PFT.max(),d.PCT_NAT_PFT.min()))
        print('land cover data for %s saved'%name)
   
if __name__=='__main__':
    create_land_cover('allgrass')
    create_land_cover('allgrass_checkerboard')
    create_land_cover('allgrass_minimal_tree')
