x0=445
x1=584
y0=1
y1=198  #lat 34.5deg
z0=1
z1=106  #z=985 m

cd /scratchu/rrolland/eNATL60-WT/2009-07

dirout='/scratchu/mclaret/eNATL60-WT_2009_07'

#fname='mesh_hgr_eNATL60SICIL_3.6'
#ncks -F -d x,$x0,$x1 -d y,$y0,$y1 ../../${fname}.nc ${dirout}/${fname}_sub.nc

#fname='ssh_T_2009-07'
#ncks -F -d x,$x0,$x1 -d y,$y0,$y1 ${fname}.nc ${dirout}/${fname}_sub.nc

run='eNATL60SICIL-BLBT02_y2009m07d'  # WT
#run='eNATL60SICIL-BLB002_y2009m08d'  # NT

#for day in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 
for day in 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do
  echo 'extracting:' ${run}${day}

  for var in votemper vosaline
  do
    fname=${run}${day}.1h_${var}
    ncks -F -d x,$x0,$x1 -d y,$y0,$y1 -d deptht,$z0,$z1 -d time_counter,1,24 ${fname}.nc ${dirout}/${fname}_sub.nc
  done

  fname=${run}${day}.1h_vozocrtx
  ncks -F -d x,$x0,$x1 -d y,$y0,$y1 -d depthu,$z0,$z1 -d time_counter,1,24 ${fname}.nc ${dirout}/${fname}_sub.nc

  fname=${run}${day}.1h_vomecrty
  ncks -F -d x,$x0,$x1 -d y,$y0,$y1 -d depthv,$z0,$z1 -d time_counter,1,24 ${fname}.nc ${dirout}/${fname}_sub.nc

  fname=${run}${day}.1h_vovecrtz
  ncks -F -d x,$x0,$x1 -d y,$y0,$y1 -d depthw,$z0,$z1 -d time_counter,1,24 ${fname}.nc ${dirout}/${fname}_sub.nc

  fname=${run}${day}.1h_sossheig
  ncks -F -d x,$x0,$x1 -d y,$y0,$y1 -d time_counter,1,24 ${fname}.nc ${dirout}/${fname}_sub.nc

done


