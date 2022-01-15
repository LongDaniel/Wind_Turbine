import numpy as np
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
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
  data[i,5]=float(value[13])*1.0

fid.close()

data0=np.genfromtxt('../acl14_clipper_2/coeff.dat')
data1=data.copy()
size0=data0.shape
print(size0)
for i in range(nline):
  data1[i,0] = (size0[0]+1+i)/250.0

data=np.concatenate((data0,data1))
print(data)

np.savetxt('coeff.dat', data)

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
plt.rc('font', family='Arial')
line_thrust =ax1.plot(data[:,0],data[:,3],label='$C_T$')
line_power =ax1.plot(data[:,0],data[:,5],label='$C_p$')
ax1.legend(bbox_to_anchor=(0.95, 0.9))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
plt.title(r'$C_T$ and $C_p$ Timehistory')
plt.xlabel('t/T')

ax1.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax1.grid(b=True, which='major', color='k', linestyle='-')
ax1.grid(b=True, which='minor', color='k', linestyle='--')
plt.savefig('fig1.png')
#plt.show()
plt.close()


