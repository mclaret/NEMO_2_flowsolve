import sys
import gsw
#import seawater as sw
import numpy as np
from netCDF4 import Dataset
from matplotlib import pyplot as plt
plt.ion()

# parameters
dirf='/scratchu/mclaret/eNATL60-WT_2009_08/'
lat=33
lon=17.5

zax_name='depth'
yax_name='yt'
xax_name='xt'

day=sys.argv[1]

fname_T="XY_d"+day+"_votemper"
fname_S="XY_d"+day+"_vosaline"

# read model output
out_T=Dataset(dirf+fname_T+'.nc',"r",format="NETCDF4")
out_S=Dataset(dirf+fname_S+'.nc',"r",format="NETCDF4")
depth=out_T.variables['depth'][:]
CT=out_T.variables['votemper'][:]
SA=out_S.variables['vosaline'][:]

# compute density
p=gsw.p_from_z(depth,lat) #TEOS10
#p=sw.pres(depth,lat)

print('computing pressure')
nt,nk,nj,ni=np.shape(SA)
p3D=np.zeros([nk,nj,ni])
for k in range(nk):
    #print(k,':',it)
    p3D[k,:,:]=p[k]

print('computing density')
# TEOS10
sigma=np.zeros([nt,nk,nj,ni])
for it in range(nt):
  sigma[it,:,:,:]=gsw.rho(SA[it,:,:,:],CT[it,:,:,:],p3D)-1000

# EOS80
#t_insitu  =np.zeros([nt,nk,nj,ni])
#rho_insitu=np.zeros([nt,nk,nj,ni])
#for it in range(nt):
#  print('time : ', it)
#  t_insitu[it,:,:,:]=sw.temp(SP[it,:,:,:],pt[it,:,:,:],p3D,0)
#  rho_insitu[it,:,:,:]=sw.dens(SP[it,:,:,:],t_insitu[it,:,:,:],p3D)

# apply mask where salinity is below 20
sigma[np.where(SA<20)]=0.e0

# write rho in-situ
out_T.close()
out_S.close()

out=Dataset(dirf+fname_T+'.nc',"a",format="NETCDF4")

sigma_nc=out.createVariable('sigma',np.float64,('time',zax_name,yax_name,xax_name), fill_value=-1.e+34)

sigma_nc[:]=sigma[:]

out.close()

