from scipy.interpolate import griddata
from netCDF4 import Dataset 
#import matplotlib.pyplot 
from matplotlib import pyplot as plt #mediterrani
import numpy as np
import xarray as xr
import gsw
import math
import sys
import os
plt.ion()

# parameters
lat=33

day=sys.argv[1]

dirf="/scratchu/mclaret/eNATL60-WT_2009_08/"
dir_ic_bc='IC_BC_v4/'

#os.mkdir(dirf+dir_ic_bc)

infile_T1=dirf+"XY_d"+day+"_votemper.nc"
infile_T2=dirf+"XY_d"+day+"_vosaline.nc"
infile_U=dirf+"XY_d"+day+"_vozocrtx.nc"
infile_V=dirf+"XY_d"+day+"_vomecrty.nc"
infile_W=dirf+"XY_d"+day+"_vovecrtz.nc"

# -------- Parameters ------

# 1st run
Lx=150.e3
Ly=Lx
ni=257; nj=257; nk=513
dz=1   # zHR grid in m
x1=211.e3 # Last  x domain location to extract
y12=130.e3  # First y domain location to extract
day0='01'

# 2nd run
#Lx=150.e3
#Ly=180.e3
#ni=257; nj=257; nk=513
#dz=1   # HR grid in m x m x m
#x1=211.e3            # Last  x domain location to extract
#y12=140.e3           # First y domain location to extract
#day0='18'

z0=-0.5     # First depth in NEMO

it1=24
# ---------------------------
# Read vars

grid_T1=Dataset(infile_T1, "r", format="NETCDF4")
grid_T2=Dataset(infile_T2, "r", format="NETCDF4")
grid_U=Dataset(infile_U, "r", format="NETCDF4")
grid_V=Dataset(infile_V, "r", format="NETCDF4")
grid_W=Dataset(infile_W, "r", format="NETCDF4")

depth_T=grid_T1.variables['depth'][:]
depth_W=grid_W.variables['depth'][:]

xxx=(grid_T1.variables['xt'][:]-grid_T1.variables['xt'][0])*1.e3  # km->m
yyy=(grid_T1.variables['yt'][:]-grid_T1.variables['yt'][0])*1.e3  # km->m
ttt=grid_T1.variables['time'][:it1]

CT=grid_T1.variables['votemper'][:it1,:,:,:]
SA=grid_T2.variables['vosaline'][:it1,:,:,:]
u=grid_U.variables['vozocrtx'][:it1,:,:,:]
v=grid_V.variables['vomecrty'][:it1,:,:,:]
w=grid_W.variables['vovecrtz'][:it1,:,:,:]

nt=np.shape(CT)[0]
# ----------------------------
# Label array dimensions

CT_xr=xr.DataArray(CT,coords={'x':xxx,'y':yyy,'z':depth_T,'t':ttt},dims=['t','z','y','x'])
SA_xr=xr.DataArray(SA,coords={'x':xxx,'y':yyy,'z':depth_T,'t':ttt},dims=['t','z','y','x'])
u_xr=xr.DataArray(u,coords={'x':xxx,'y':yyy,'z':depth_T,'t':ttt},dims=['t','z','y','x'])
v_xr=xr.DataArray(v,coords={'x':xxx,'y':yyy,'z':depth_T,'t':ttt},dims=['t','z','y','x'])
w_xr=xr.DataArray(w,coords={'x':xxx,'y':yyy,'z':depth_W,'t':ttt},dims=['t','z','y','x'])

del CT,SA,u,v,w
# ----------------------------
# Extract domain of interest and Interpolate into HR coordinates at the same time

Lz=dz*(nk-1)*(-1)-dz  # First flow_solve depth is z=-dz in NEMO
#Lx=ni*dx
#Ly=nj*dy

xHR=np.linspace(x1-Lx,x1,num=ni)
yHR=np.linspace(y12-Ly/2,y12+Ly/2,num=nj)
zHR=np.linspace(Lz,-dz,num=nk) 

# Convert coordinates to flow_solve convention (names and zero start)
x_fs=xHR-np.min(xHR)
y_fs=yHR-np.min(yHR)
z_fs=zHR-np.min(zHR)  # z=0 at k=1 like in atmos

# ----------------------------
# Depth in 3D
# In the TEOS10 used in NEMO, pressure is approcimate by depth in meters
z3D=np.zeros([nk,nj,ni])
for k in range(nk):
    z3D[k,:,:]=zHR[k]*(-1)

# ----------------------------

time_attrs = dict(standard_name="time",
                  units="seconds since 1900-01-01 00:00:00",
                  calendar="gregorian")

u_HR=np.zeros((1,nk,nj,ni))
v_HR=np.zeros((1,nk,nj,ni))
w_HR=np.zeros((1,nk,nj,ni))


for it in range(nt):
  print('regridding to flow_solve grid time T/S/u/v/w. dd:hh', day,':',it)
  CT_HR=CT_xr.isel(t=it).interp(x=xHR,y=yHR,z=zHR,method='linear').values
  SA_HR=SA_xr.isel(t=it).interp(x=xHR,y=yHR,z=zHR,method='linear').values
  u_HR[0,:,:,:]=u_xr.isel(t=it).interp(x=xHR,y=yHR,z=zHR,method='linear').values
  v_HR[0,:,:,:]=v_xr.isel(t=it).interp(x=xHR,y=yHR,z=zHR,method='linear').values
  w_HR[0,:,:,:]=w_xr.isel(t=it).interp(x=xHR,y=yHR,z=zHR,method='linear').values

  sig_HR=gsw.rho(SA_HR,CT_HR,z3D)-1000
  del SA_HR, CT_HR

  sig_anom=np.zeros([1,nk,nj,ni])
  if (day==day0):
    if it==0:
      sig_mean=np.mean(np.mean(sig_HR,axis=2),axis=1)

    for k in range(nk):
      sig_anom[0,k,:,:]=sig_HR[k,:,:]-sig_mean[k]
  else:
    if it==0:
      data=Dataset(dirf+dir_ic_bc+'XYZ_000.nc','r',format="NETCDF4")
      sig_mean=data.variables['s1_bar'][:]-1000

    for k in range(nk):
      sig_anom[0,k,:,:]=sig_HR[k,:,:]-sig_mean[k]

  ds = xr.Dataset(
  data_vars=dict(
    s1_bar=(["z"], sig_mean+1000, {'units':'kg/m3', 'missing_value': -1.e+34}),
    s1 =(["time","z","y","x"], sig_anom,{'units':'kg/m3', 'missing_value': -1.e+34}),
    u  =(["time","z","y","x"], u_HR ,{'units':'m/s', 'missing_value': -1.e+34}),
    v  =(["time","z","y","x"], v_HR ,{'units':'m/s', 'missing_value': -1.e+34}),
    w  =(["time","z","y","x"], w_HR ,{'units':'m/s', 'missing_value': -1.e+34}),
  ),
  coords=dict(
    x=(["x"], x_fs),
    y=(["y"], y_fs),
    z=(["z"], z_fs),
    time=(["time"],ttt[it:it+1],time_attrs),
  ),
  attrs=dict(description="NEMO_2_FLOWSOLVE"),)

  if day==day0 and it==0:
    # Create files
    ds.to_netcdf(dirf+dir_ic_bc+'XYZ_000.nc')   # Initial condition file

  if it<10:
    itstr='0'+str(it)
  else:
    itstr=str(it)

  # Boundary files
  ds.isel(x=0).to_netcdf(dirf+dir_ic_bc+'west_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')
  ds.isel(x=ni-1).to_netcdf(dirf+dir_ic_bc+'east_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')

  ds.isel(y=0).to_netcdf(dirf+dir_ic_bc+'south_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')
  ds.isel(y=nj-1).to_netcdf(dirf+dir_ic_bc+'north_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')

  ds.isel(z=0).to_netcdf(dirf+dir_ic_bc+'bottom_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')
  ds.isel(z=nk-1).to_netcdf(dirf+dir_ic_bc+'top_d'+day+'_it'+itstr+'.nc',unlimited_dims='time')
