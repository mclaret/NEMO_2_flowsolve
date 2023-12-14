#MClaret@LOCEAN@Sep'23

from scipy.interpolate import griddata
from netCDF4 import Dataset
#import matplotlib.pyplot 
from matplotlib import pyplot as plt #mediterrani
import numpy as np
import math
plt.ion()

#################

print('reading files')

# Time indices of the regridding
it0=0; it1=408
dirf="/scratchu/mclaret/eNATL60-WT_2009_07/"

grid_T=Dataset(dirf+"ssh_T_2009-07_sub.nc", "r", format="NETCDF4")
coord=Dataset("/scratchu/mclaret/eNATL60-WT_2009_08/"+"mesh_hgr_eNATL60SICIL_3.6_sub.nc","r", format="NETCDF4")

print(grid_T.dimensions.keys())
print(grid_T.variables.keys())

# read vars
lat_T=grid_T.variables['nav_lat'][:] #NOTE: a colon is needed otherwise is shows how the data is
lon_T=grid_T.variables['nav_lon'][:]

delx_T=np.squeeze(coord.variables['e1t'][:])
delx_U=np.squeeze(coord.variables['e1u'][:])
delx_V=np.squeeze(coord.variables['e1v'][:])

dely_T=np.squeeze(coord.variables['e2t'][:])
dely_U=np.squeeze(coord.variables['e2u'][:])
dely_V=np.squeeze(coord.variables['e2v'][:])

time_count=grid_T.variables['time_counter'][it0:it1]

ssh=grid_T.variables['sossheig'][it0:it1,:,:]

grid_T.close()
coord.close()

nt,ny,nx=np.shape(ssh)

print(it0,':',it1,':',nt)


#####################################################3
## Coordinates switch Geographical -> Cartesian

#--Uncomment if coordinates file is not available--#

#pi=4.*math.atan(1.)
#deglen=111198.6     #  SCALAR one degree latitude in meters

## Tgrid and Ugrid
#del_lat_T=np.diff(lat_T,n=1,axis=0)
#deglenx_T=deglen*np.cos(pi*lat_T/180)       #2D one deg lon at latitude y 
#delx_T=deglenx_T*1/64 
#dely_T=deglen*del_lat_T

#xx_T=np.zeros((ny,nx))
#xx =np.cumsum(delx_T[:,0:376],axis=-1)/1.e3  # Tgrid in km
#xx_T[:,1:377]=xx

#yy_T=np.zeros((ny,nx))
#yy=np.cumsum(dely_T,axis=0)/1.e3   # new Y grid in km
#yy_T[1:170,:]=yy

#----------------------------------------------------

print('Build Cartesian grids')

xx_T=np.zeros((ny,nx))
xx =np.cumsum(delx_U[:,0:-1],axis=-1)/1.e3  # Tgrid in km
xx_T[:,1:nx]=xx

yy_T=np.zeros((ny,nx))
yy=np.cumsum(dely_V[0:-1,:],axis=0)/1.e3   # new Y grid in km
yy_T[1:ny,:]=yy

del_regrid=1.5  # regridded into 1.5km evenly spaced grid along x and y
nx2=int(np.max(xx[ny-1,:])/del_regrid)  #np.max(xx[ny-1,:])
ny2=int(np.max(yy[:,nx-1])/del_regrid)  #np.max(yy[:,nx-1])

x2v=np.linspace(0,nx2,nx2+1)*del_regrid
y2v=np.linspace(0,ny2,ny2+1)*del_regrid
xx2,yy2=np.meshgrid(x2v,y2v)

xx_T_vect=np.concatenate(np.transpose(xx_T))
yy_T_vect=np.concatenate(np.transpose(yy_T))

#####################################################3

print('prepare output file and write coordinates')
outxy=Dataset(dirf+"XY_m07.ssh.nc",'w',format='NETCDF4')
# create dimensions
xdim = outxy.createDimension('xt', nx2+1)
ydim = outxy.createDimension('yt', ny2+1)
tdim = outxy.createDimension('time', None)
print('dimension variables:')
for dimname in outxy.dimensions.keys():     
    dim = outxy.dimensions[dimname]    
    print ( dimname, len(dim), dim.isunlimited() )
# create variables
x_nc=outxy.createVariable('xt',np.float32,('xt'), fill_value=-1.e+34)
y_nc=outxy.createVariable('yt',np.float32,('yt'), fill_value=-1.e+34)
t_nc=outxy.createVariable('time',np.float64,('time'))
t_nc.units="seconds since 1900-01-01 00:00:00"
t_nc.time_origin="1900-01-01 00:00:00"
# write dimensions
x_nc[:]=x2v[:]
y_nc[:]=y2v[:]
t_nc[:]=time_count[:]
# create variables
ssh_nc = outxy.createVariable('SSH', np.float32, ('time','yt','xt'), fill_value=-1.e+34)
print('variables:')
for varname in outxy.variables.keys():    
    var_nc = outxy.variables[varname]    
    print (varname, var_nc.dtype, var_nc.dimensions, var_nc.shape)

#####################################################

print('regridding:', varname)

# REGRID SSH
for t in range(nt):
  print('regridding T times,', t)
  tmp1=np.concatenate(np.transpose(ssh[t,:,:]))
  tmp2=griddata((xx_T_vect,yy_T_vect),tmp1,(np.transpose(xx2),np.transpose(yy2)))
  tmp3=np.transpose(tmp2)
  ssh_nc[t,:,:]=tmp3

# close writing file
outxy.close()
