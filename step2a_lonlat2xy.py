#MClaret@LOCEAN@Sep'23. 
#Python 3.7.4 (default, Sep  7 2019, 19:52:29) 
#port select --list python
#sudo port select --set python python37

from scipy.interpolate import griddata
from netCDF4 import Dataset
#import matplotlib.pyplot 
from matplotlib import pyplot as plt #mediterrani
import numpy as np
import math
import sys
plt.ion()

#################

print('reading files')

# Time indices of the regridding
it0=0; 
#it1=24;
dirf="/scratchu/mclaret/eNATL60-WT_2009_07/"

varname=sys.argv[1]
day=sys.argv[2]

if (varname=='votemper') or (varname=='vosaline'):
  zname='deptht'
elif (varname=='vozocrtx'):
  zname='depthu'
elif (varname=='vomecrty'):
  zname='depthv'
elif (varname=='vovecrtz'):
  zname='depthw'

infile="eNATL60SICIL-BLBT02_y2009m07d"+day+".1h_"+varname+"_sub.nc"
grid=Dataset(dirf+infile, "r", format="NETCDF4")
coordH=Dataset("/scratchu/mclaret/eNATL60-WT_2009_08/"+"mesh_hgr_eNATL60SICIL_3.6_sub.nc","r", format="NETCDF4")

#print(grid.dimensions.keys())
#print(grid.variables.keys())

# read vars
print('reading: ', infile)
lat=grid.variables['nav_lat'][:] #NOTE: a colon is needed otherwise is shows how the data is
lon=grid.variables['nav_lon'][:]

delx_T=np.squeeze(coordH.variables['e1t'][:])
delx_U=np.squeeze(coordH.variables['e1u'][:])
delx_V=np.squeeze(coordH.variables['e1v'][:])

dely_T=np.squeeze(coordH.variables['e2t'][:])
dely_U=np.squeeze(coordH.variables['e2u'][:])
dely_V=np.squeeze(coordH.variables['e2v'][:])

time_count=grid.variables['time_counter'][:]

var=grid.variables[varname][:,:,:,:]
depth=np.squeeze(grid.variables[zname][:])

grid.close()
coordH.close()

nt,nz,ny,nx=np.shape(var)

it1=nt
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

## Vgrid
#lat_V2=np.zeros((ny+1,nx))
#lat_V2[0,:]=lat_T[0,:]
#lat_V2[1:171,:]=lat_V
#del_lat_V=np.diff(lat_V2,n=1,axis=0)
#deglenx_V=deglen*np.cos(pi*lat_V/180)       #2D one deg lon at latitude y 
#delx_V=deglenx_V*1/64 
#dely_V=deglen*del_lat_V

#xx_T=np.zeros((ny,nx))
#xx_V=np.zeros((ny,nx))
#xx =np.cumsum(delx_T[:,0:376],axis=-1)/1.e3  # Tgrid in km
#xx2=np.cumsum(delx_V[:,0:376],axis=-1)/1.e3  # Tgrid in km
#xx_T[:,1:377]=xx
#xx_U=xx_T+delx_T/2/1.e3              # Ugrid in km  Ugrid is displaced delx/2 relative to Tgrid
#xx_V[:,1:377]=xx2

#yy_T=np.zeros((ny,nx))
#yy=np.cumsum(dely_T,axis=0)/1.e3   # new Y grid in km
#yy_T[1:170,:]=yy
#yy_U=yy_T
#yy_V=np.cumsum(dely_V,axis=0)/1.e3   # new Y grid in km

#--------------------------------------------------------
# BUILD GRIDS

print('Build Cartesian grids')

# Build NEMO x-grids
xx_T=np.zeros((ny,nx))
xx_U=np.zeros((ny,nx))
xx_V=np.zeros((ny,nx))

xx0=np.cumsum(delx_U[:,0:-1],axis=-1)/1.e3   # Tgrid in km
xx1=np.cumsum(delx_T[:,1:nx],axis=-1)/1.e3   # Ugrid in km
delx_V2=0.5*(delx_V[:,0:nx-1]+delx_V[:,1:nx]) # delx at U-points at V-points
xx2=np.cumsum(delx_V2,axis=-1)/1.e3   # Vgrid in km

xx_T[:,1:nx]=xx0
xx_U[:,1:nx]=xx1
for j in range(0,ny):
  xx_U[j,:]=xx_U[j,:]+delx_T[j,0]/2/1.e3  # Ugrid is displaced delx/2 relative to Tgrid

xx_V[:,1:nx]=xx2

# Build NEMO y-grids
yy_T=np.zeros((ny,nx))
yy_U=np.zeros((ny,nx))
yy_V=np.zeros((ny,nx))

yy0=np.cumsum(dely_V[0:-1],axis=0)/1.e3       # new Y grid in km
dely_V2=0.5*(dely_U[1:ny,:]+dely_U[1:ny,:]) # delx at U-points at V-points
yy1=np.cumsum(dely_V2,axis=0)/1.e3
yy2=np.cumsum(dely_T[1:ny,:],axis=0)/1.e3

yy_T[1:ny,:]=yy0
yy_U[1:ny,:]=yy1
yy_V[1:ny,:]=yy2   
for i in range(0,nx):
  yy_V[:,i]=yy_V[:,i]+dely_T[0,i]/2/1.e3 # Vgrid is displaced dely/2 relative to Tgrid

# Build Cartesian grid
del_regrid=1.5  # regridded into 1.5km evenly spaced grid along x and y
nx2=int(np.max(xx_T[ny-1,:])/del_regrid)  #np.max(xx[ny-1,:])
ny2=int(np.max(yy_T[:,nx-1])/del_regrid)  #np.max(yy[:,nx-1])

x2v=np.linspace(0,nx2,nx2+1)*del_regrid
y2v=np.linspace(0,ny2,ny2+1)*del_regrid
xx2,yy2=np.meshgrid(x2v,y2v)

xx_T_vect=np.concatenate(np.transpose(xx_T))
yy_T_vect=np.concatenate(np.transpose(yy_T))
xx_U_vect=np.concatenate(np.transpose(xx_U))
yy_U_vect=np.concatenate(np.transpose(yy_U))
xx_V_vect=np.concatenate(np.transpose(xx_V))
yy_V_vect=np.concatenate(np.transpose(yy_V))

#####################################################3

outfile=dirf+"XY_d"+day+"_"+varname+".nc"
print('prepare output file and write coordinates: ', outfile)
outxy=Dataset(outfile,'w',format='NETCDF4')
# create dimensions
xdim = outxy.createDimension('xt', nx2+1)
ydim = outxy.createDimension('yt', ny2+1)
zdim = outxy.createDimension('depth', nz)
tdim = outxy.createDimension('time', None)
print('dimension variables:')
for dimname in outxy.dimensions.keys():     
    dim = outxy.dimensions[dimname]    
    print ( dimname, len(dim), dim.isunlimited() )

# create variables
x_nc=outxy.createVariable('xt',np.float32,('xt'), fill_value=-1.e+34)
y_nc=outxy.createVariable('yt',np.float32,('yt'), fill_value=-1.e+34)
z_nc=outxy.createVariable('depth',np.float32,('depth'), fill_value=-1.e+34)
t_nc=outxy.createVariable('time',np.float64,('time'))
z_nc.positive="up"
t_nc.units="seconds since 1900-01-01 00:00:00"
t_nc.time_origin="1900-01-01 00:00:00"
# write dimensions
x_nc[:]=x2v[:]
y_nc[:]=y2v[:]
z_nc[:]=np.flip(depth[:])*(-1)
#z_nc[:]=np.linspace(0,nz,nz)+1
t_nc[:]=time_count[:]
# create variables
var_nc = outxy.createVariable(varname, np.float32, ('time','depth','yt','xt'), fill_value=-999)
print (varname, var_nc.dtype, var_nc.dimensions, var_nc.shape)

#####################################################

print('regridding: ', varname)

# REGRID from T-grid
if (varname=='votemper') or (varname=='vosaline') or (varname=='vovecrtz') or (varname=='votkeavt'):
  for t in range(nt):
    print('regridding T/S/W/Kz times,', t)
    for k in range(nz):
        tmp1=np.concatenate(np.transpose(var[t,k,:,:]))
        tmp2=griddata((xx_T_vect,yy_T_vect),tmp1,(np.transpose(xx2),np.transpose(yy2)))
        tmp3=np.transpose(tmp2)
        var_nc[t,nz-1-k,:,:]=tmp3   # flip vertical levels

# REGRID from U-grid
if (varname=='vozocrtx'):
  for t in range(nt):
    print('regridding U times,', t)
    for k in range(nz):
      tmp1=np.concatenate(np.transpose(var[t,k,:,:]))
      tmp2=griddata((xx_U_vect,yy_U_vect),tmp1,(np.transpose(xx2),np.transpose(yy2)))
      tmp3=np.transpose(tmp2)
      var_nc[t,nz-1-k,:,:]=tmp3    # flip vertical levels

# REGRID from V-velocity
if (varname=='vomecrty'):
  for t in range(nt):
    print('regridding V times,', t)
    for k in range(nz):
      tmp1=np.concatenate(np.transpose(var[t,k,:,:]))
      tmp2=griddata((xx_V_vect,yy_V_vect),tmp1,(np.transpose(xx2),np.transpose(yy2)))
      tmp3=np.transpose(tmp2)
      var_nc[t,nz-1-k,:,:]=tmp3    # flip vertical levels

# close writing file
outxy.close()
