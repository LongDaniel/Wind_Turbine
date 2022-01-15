import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import tecplot_io as tec

#os.makedirs("stats", exist_ok = True)
os.system("mkdir -p stats")

## motion_ut1
config = 2 
icase_s = 128 
icase_e = 160 + 1

## emptyhos
#config = 3 
#icase_s = 11 
#icase_e = 16 + 1

## empty_r1
#config = 4 
#icase_s = 81 
#icase_e = 99 + 1

if config == 2:
  nvar = 8
  nz = 65 
  f_pre = "motion_ut2"
  icopy_basic = True
  icopy_wtforce = True 
  icopy_wakeline = True

ncase = icase_e - icase_s

flux_stat = np.zeros(ncase)

for ii in range(ncase):

  icase = icase_s + ii
  fname = "./"+f_pre+"_" + str(icase)+"/fort.930.dat"
  print("Reading "+fname)
  
  data0 = tec.tecplot_reader(fname, [nz, 1, 1, nvar], 2)
  data1 = data0[:,0,0,:].reshape(nz, nvar)
  
  flux_stat[ii] = np.trapz(data1[:,1], data1[:,0])
  
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(1,1,1)
  ax1.plot(data1[:,1], data1[:,0])
  plt.xlim([0, 1.4])
  plt.ylim([0, 0.46])
  plt.savefig("um_{0:03d}.png".format(icase))
  plt.close()
  
  # copy err1, fort.930.dat
  if icopy_basic:
    tmp_dir = "./stats/"+f_pre+"_"+str(icase)
    #os.makedirs(tmp_dir, exist_ok = True)
    os.system("mkdir -p "+tmp_dir)
    tmp_cmd = "cp -t "+tmp_dir+" "+f_pre+"_"+str(icase)+"/err1"+" "+f_pre+"_"+str(icase)+"/fort.930.dat"
    print("Executing: "+tmp_cmd)
    os.system(tmp_cmd)
  
  # copy Turbine_AL_001.dat (001~016), Nacelle_001.dat (001~048) 
  if icopy_wtforce:
    tmp_dir = "./stats/"+f_pre+"_"+str(icase)
    #os.makedirs(tmp_dir, exist_ok = True)
    os.system("mkdir -p "+tmp_dir)
    tmp_cmd = "cp -t "+tmp_dir+" "+f_pre+"_"+str(icase)+"/Turbine_AL_*.dat"+" "+f_pre+"_"+str(icase)+"/Nacelle_*.dat"
    print("Executing: "+tmp_cmd)
    os.system(tmp_cmd)
  
  # copy wakeline
  if icopy_wakeline:
    is_wake = False
    try:
      with open(f_pre+"_"+str(icase)+"/wakeline_001.dat") as file:
        is_wake = True
        tmp_cmd = "cp -t "+tmp_dir+" "+f_pre+"_"+str(icase)+"/wakeline_*"
        print("Executing: "+tmp_cmd)
        os.system(tmp_cmd)
    except IOError:
      print("No wakeline log in directory: "+f_pre+"_"+str(icase))

fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)
ax2.plot(flux_stat, '+')
plt.savefig("flux_stat.png")
plt.close()

print("\n----------Flux stats:-------------\nicase\tflux")
for ii in range(ncase):
  icase = icase_s + ii
  print(str(icase)+"\t"+str(flux_stat[ii]))

