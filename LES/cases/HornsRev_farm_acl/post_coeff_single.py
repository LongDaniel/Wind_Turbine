import numpy as np
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

import os

foldername = os.path.relpath(".","..")
print(foldername)

dt = 0.189138683933 
nturbinex = 3
nturbiney = 3
nturbine = nturbinex * nturbiney

## import thrust coefficient log
fid = open('log.cthrust','r')
content=fid.readlines()
# Thrust=  0.10233803188671066      , Torque=   6.1668522808625047E-003 , C_Thrust=  0.38109671967396186      , Power=   6.3533464098887757E-002 , C_Power=  0.12840488196034763     

nline = len(content)
nt = nline / nturbine
data1=np.zeros((nturbine*nt,6))

for i in range(nt):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data1[i2,0]=(i+1)*dt
    data1[i2,1]=value[1] # thrust
    data1[i2,2]=value[4] # torque
    data1[i2,3]=value[7] # c_thrust
    data1[i2,4]=value[10] # power
    data1[i2,5]=float(value[13])*1.0 # c_power

fid.close()


np.savetxt('coeff.dat', data1)


## import reference velocity log
fid = open('log.uref','r')
content=fid.readlines()
# Turbine_           1 :angvel=   50.735637598099999      , TSR=   4.4059452767309510      , Uref=  0.86364504796591524

nline = len(content)
data2=np.zeros((nline,4))

for i in range(nt):
  for j in range(nturbine):
    i2 = nturbine * i + j
    value=content[i2].split()
    #print(value)
    data2[i2,0]=(i+1)*dt
    data2[i2,1]=value[3] # angvel
    data2[i2,2]=value[6] # TSR
    data2[i2,3]=value[9] # Uref

fid.close()

np.savetxt('uref.dat', data2)


