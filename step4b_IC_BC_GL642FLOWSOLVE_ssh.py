from scipy.interpolate import griddata
from netCDF4 import Dataset 
#import matplotlib.pyplot 
from matplotlib import pyplot as plt #mediterrani
import numpy as np
import xarray as xr
import math
import sys
plt.ion()

day=sys.argv[1]

dirf="/scratchu/mclaret/eNATL60-WT_2009_08/"
dir_ic_bc='IC_BC_v4/'
day0='01'
#day0='18'

#infile=dirf+"XY_m07.ssh_18_31.nc"
infile=dirf+"XY_m08.ssh.nc"

iday=int(day)-1-(int(day0)-1)
it0=iday*24
it1=it0+24

# -------- Parameters ------

# First run
Lx=150.e3
Ly=Lx
ni=257; nj=257
x1=211.e3            # Last  x domain location to extract
y12=130.e3           # First y domain location to extract

# Second run
#Lx=150.e3
#Ly=180.e3
#ni=257; nj=257
#x1=211.e3            # Last  x domain location to extract
#y12=140.e3           # First y domain location to extract

z0=-0.5              # First depth in NEMO

# ---------------------------
# Read vars

grid_T=Dataset(infile, "r", format="NETCDF4")

xxx=(grid_T.variables['xt'][:]-grid_T.variables['xt'][0])*1.e3  # km->m
yyy=(grid_T.variables['yt'][:]-grid_T.variables['yt'][0])*1.e3  # km->m
ttt=grid_T.variables['time'][it0:it1]

ssh=grid_T.variables['SSH'][it0:it1,:,:]

nt=np.shape(ssh)[0]
# ----------------------------
# Label array dimensions

ssh_xr=xr.DataArray(ssh,coords={'x':xxx,'y':yyy,'t':ttt},dims=['t','y','x'])

# ----------------------------
# Extract domain of interest and Interpolate into HR coordinates at the same time

#Lx=ni*dx
#Ly=nj*dy

print('regridding ssh day', day)

xHR=np.linspace(x1-Lx,x1,num=ni)
yHR=np.linspace(y12-Ly/2,y12+Ly/2,num=nj)

ssh_HR_xr=ssh_xr.interp(x=xHR,y=yHR,method='linear')

# ----------------------------
# Convert coordinates to flow_solve convention (names and zero start)
x_fs=xHR-np.min(xHR)
y_fs=yHR-np.min(yHR)

time_attrs = dict(standard_name="time",
                  units="seconds since 1900-01-01 00:00:00",
                  calendar="gregorian")

ds = xr.Dataset(
data_vars=dict(
    eta =(["time","y","x"], ssh_HR_xr.values,{'units':'m', 'missing_value': -1.e+34}),
),
coords=dict(
    x=(["x"], x_fs),
    y=(["y"], y_fs),
    time=(["time"],ttt,time_attrs),
),
attrs=dict(description="NEMO_2_FLOWSOLVE"),)

# Boundary condition files
ds.to_netcdf(dirf+dir_ic_bc+'top_d'+day+'.nc', mode="a")



