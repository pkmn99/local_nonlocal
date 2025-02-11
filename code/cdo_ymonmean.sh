# use cdo to compute multi-year monthly mean
path=../data/result/
cdo ymonmean ${path}ibm.TSc.local_nonlocal.mon.nc ${path}ibm.TSc.local_nonlocal.ymonmean.nc
cdo ymonmean ${path}offline.TSc.local_nonlocal.mon.nc ${path}offline.TSc.local_nonlocal.ymonmean.nc
cdo ymonmean ${path}Tair.TSc.local_nonlocal.mon.nc ${path}Tair.TSc.local_nonlocal.ymonmean.nc
cdo ymonmean ${path}interpolation.TSc.local_nonlocal.mon.nc ${path}interpolation.TSc.local_nonlocal.ymonmean.nc
cdo ymonmean ${path}subgrid_offline.TSc.local_nonlocal.mon.nc ${path}subgrid_offline.TSc.local_nonlocal.ymonmean.nc
cdo ymonmean ${path}subgrid.TSc.local_nonlocal.mon.nc ${path}subgrid.TSc.local_nonlocal.ymonmean.nc
echo "done"
