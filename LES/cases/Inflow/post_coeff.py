import numpy as np
import matplotlib
# Force matplotlib to not use any Xwindows backend
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

fid = open('log.cthrust','r')
content=fid.readlines()
# Thrust=  0.10233803188671066      , Torque=   6.1668522808625047E-003 , C_Thrust=  0.38109671967396186      , Power=   6.3533464098887757E-002 , C_Power=  0.12840488196034763     

nline = len(content)
data=np.zeros((nline,6))

for i in range(nline):
  value=content[i].split()
  #print(value)
  data[i,0]=(i+1)/250.0
  data[i,1]=value[1]
  data[i,2]=value[4]
  data[i,3]=value[7]
  data[i,4]=value[10]
  data[i,5]=value[13]

np.savetxt('coeff.dat', data)

#plt.rc('text', usetex=True)
plt.rc('font', family='serif')
line_thrust, =plt.plot(data[:,0],data[:,3],label='$C_T$')
line_power, =plt.plot(data[:,0],data[:,5],label='$C_p$')
plt.legend(handles=[line_thrust,line_power])
plt.axis([0,40,0,1])
plt.title(r'$C_T$ and $C_p$ Timehistory')
plt.xlabel('t/T')
plt.savefig('fig1.png')
plt.show()
plt.close()


