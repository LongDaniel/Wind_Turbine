import numpy as np
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

import os

foldername = os.path.relpath(".","..")
print(foldername)

dt = 0.001

## import thrust coefficient log
fid = open('log.ud','r')
content=fid.readlines()
#          1 ud_instant=  0.59557850926407163      , ud_mean=  0.60611835966224104

nline = len(content)
data1=np.zeros((nline,3))

for i in range(nline):
  value=content[i].split()
  #print(value)
  data1[i,0]=(i+1)*dt
  data1[i,1]=value[2] # ud_instant 
  data1[i,2]=value[5] # ud_mean

fid.close()

'''
i_win1 = 2000
t_win1 = i_win1 * dt 
mean1 = data1.copy()
for i in range(nline):
  mean1[i,0] = data1[i,0]
  if i<i_win1:
    for j in range(5):
      mean1[i,j+1] = np.average(data1[0:(i+1),j+1])
  else:
    mean1[i,1:6] = data1[i,1:6]*(1.0/i_win1)+mean1[i-1,1:6]*(1.0-1.0/i_win1)
'''

#data0=np.genfromtxt('../acl14_clipper_2/coeff.dat')
#data1=data.copy()
#size0=data0.shape
#print(size0)
#for i in range(nline):
#  data1[i,0] = (size0[0]+1+i)/250.0
#
#data=np.concatenate((data0,data1))
#print(data)

#np.savetxt('coeff.dat', data1)
np.savetxt('uref.dat', data1)

'''
## import reference velocity log
fid = open('log.uref','r')
content=fid.readlines()
# Turbine_           1 :angvel=   50.735637598099999      , TSR=   4.4059452767309510      , Uref=  0.86364504796591524

nline = len(content)
data2=np.zeros((nline,4))

for i in range(nline):
  value=content[i].split()
  #print(value)
  data2[i,0]=(i+1)*0.001/0.123841654597
  data2[i,1]=value[3] # angvel
  data2[i,2]=value[6] # TSR
  data2[i,3]=value[9] # Uref

fid.close()

i_win2 = 1000
t_win2 = i_win2 * dt 
mean2 = data2.copy()
for i in range(nline):
  mean2[i,0] = data2[i,0]
  if i<i_win2:
    for j in range(3):
      mean2[i,j+1] = np.average(data2[0:(i+1),j+1])
  else:
    mean2[i,1:4] = data2[i,1:4]*(1.0/i_win2)+mean2[i-1,1:4]*(1.0-1.0/i_win2)
'''

#np.savetxt('uref.dat', data1)

'''
## plot C_Thrust and C_Power
fig1 = plt.figure()
ax1_1 = fig1.add_subplot(1,1,1)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', family='Arial')
line_thrust =ax1_1.plot(data1[:,0],data1[:,3],':',label='$C_T$')
line_thrustmean = ax1_1.plot(mean1[:,0], mean1[:,3],'-',label=r'$\left<C_T\right>_{t}$')
line_power =ax1_1.plot(data1[:,0],data1[:,5],'y-.',label='$C_p$')
line_powermean = ax1_1.plot(mean1[:,0], mean1[:,5],'--',label=r'$\left<C_p\right>_{t}$')
ax1_1.legend(bbox_to_anchor=(0.95, 0.9))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$C_T$ and $C_p$ time-history')
plt.xlabel('$t/T$')
plt.ylabel('$C_T,\ C_p$')

ax1_1.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax1_1.grid(b=True, which='major', color='k', linestyle='-')
ax1_1.grid(b=True, which='minor', color='k', linestyle='--')
plt.savefig('fig1_ct_cp_'+foldername+'.png')
plt.savefig('fig1_ct_cp_'+foldername+'.pdf', format='pdf')
#plt.show()
plt.close()
'''

## plot Uref 
fig2 = plt.figure()
ax2_1 = fig2.add_subplot(1,1,1)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', family='Arial')
line_uref = ax2_1.plot(data1[:,0],data1[:,1],'--',label='$U_{ref}$')
line_urefmean = ax2_1.plot(data1[:,0], data1[:,2],'-',label=r'$\left<U_{ref}\right>_{t}$')
ax2_1.legend(bbox_to_anchor=(0.95, 0.9))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$U_{ref}$ time-history')
#plt.xlabel('t/T')
plt.ylabel('$U_{ref}$')

ax2_1.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2_1.grid(b=True, which='major', color='k', linestyle='-')
ax2_1.grid(b=True, which='minor', color='k', linestyle='--')

'''
ax2_2 = fig2.add_subplot(2,1,2)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', family='Arial')
line_TSR = ax2_2.plot(data2[:,0],data2[:,2],'--',label='$TSR$')
line_TSRmean = ax2_2.plot(mean2[:,0], mean2[:,2],'-',label=r'$\left<TSR\right>_{t}$')
ax2_2.legend(bbox_to_anchor=(0.95, 0.9))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$TSR$ time-history')
plt.xlabel('$t/T$')
plt.ylabel('$TSR$')

ax2_2.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2_2.grid(b=True, which='major', color='k', linestyle='-')
ax2_2.grid(b=True, which='minor', color='k', linestyle='--')
'''

plt.savefig('fig2_uref_tsr_'+foldername+'.png')
plt.savefig('fig2_uref_tsr_'+foldername+'.pdf', format='pdf')
#plt.show()
plt.close()

