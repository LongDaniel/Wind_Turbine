import numpy as np
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

import os

foldername = os.path.relpath(".","..")
print(foldername)

'''
  data1[i,0]=(i+1)*0.001/0.123841654597/2
  data1[i,1]=value[1] # thrust
  data1[i,2]=value[4] # torque
  data1[i,3]=value[7] # c_thrust
  data1[i,4]=value[10] # power
  data1[i,5]=float(value[13])*1.0 # c_power
'''

ncase = 4 
dt = 0.00136225820057 

temp0 = np.array([])
for i in range(ncase):
  temp0 = temp0.copy()
  temp1 = np.genfromtxt('coeff_'+str(i+1)+'.dat')
  if i == 0:
    temp0 = temp1.copy()
  else:
    temp0 = np.concatenate((temp0,temp1))

data1 = temp0.copy()
size1 = temp0.shape
for i in range(size1[0]):
  #data1[i,0] = (i+1)*dt/0.123841654597/2.0
  data1[i,0] = i*dt

temp0 = np.zeros(0)
for i in range(ncase):
  temp0 = temp0.copy()
  temp1 = np.genfromtxt('uref_'+str(i+1)+'.dat')
  if i==0:
    temp0 = temp1.copy()
  else:
    temp0 = np.concatenate((temp0,temp1))

data2 = temp0.copy()
size2 = temp0.shape
for i in range(size2[0]):
  #data2[i,0] = (i+1)*dt/0.123841654597
  data2[i,0] = i*dt

#print(data)

np.savetxt('coeff_all.dat', data1)
np.savetxt('uref_all.dat', data2)

i_win1 = 2000
t_win1 = i_win1 * dt 
mean1 = data1.copy()
nline = (data1.shape)[0]
for i in range(nline):
  mean1[i,0] = data1[i,0]
  if i<i_win1:
    for j in range(5):
      mean1[i,j+1] = np.average(data1[0:(i+1),j+1])
  else:
    mean1[i,1:6] = data1[i,1:6]*(1.0/i_win1)+mean1[i-1,1:6]*(1.0-1.0/i_win1)

i_win2 = 2000
t_win2 = i_win2 * dt 
mean2 = data2.copy()
nline = (data2.shape)[0]
for i in range(nline):
  mean2[i,0] = data2[i,0]
  if i<i_win2:
    for j in range(3):
      mean2[i,j+1] = np.average(data2[0:(i+1),j+1])
  else:
    mean2[i,1:4] = data2[i,1:4]*(1.0/i_win2)+mean2[i-1,1:4]*(1.0-1.0/i_win2)


'''
  data2[i,0]=(i+1)*0.001/0.123841654597
  data2[i,1]=value[3] # angvel
  data2[i,2]=value[6] # TSR
  data2[i,3]=value[9] # Uref
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
ax1_1.legend(bbox_to_anchor=(0.2, 0.95))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$C_T$ and $C_p$ time-history')
plt.xlabel('$t$')
plt.ylabel('$C_T,\ C_p$')

ax1_1.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax1_1.grid(b=True, which='major', color='k', linestyle='-')
ax1_1.grid(b=True, which='minor', color='k', linestyle='--')
plt.savefig('fig1_ct_cp_'+foldername+'.png')
plt.savefig('fig1_ct_cp_'+foldername+'.pdf', format='pdf')
#plt.show()
plt.close()


## plot Uref and TSR
fig2 = plt.figure()
ax2_1 = fig2.add_subplot(2,1,1)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', family='Arial')
line_uref = ax2_1.plot(data2[:,0],data2[:,3],'--',label='$U_{ref}$')
line_urefmean = ax2_1.plot(mean2[:,0], mean2[:,3],'-',label=r'$\left<U_{ref}\right>_{t}$')
ax2_1.legend(bbox_to_anchor=(0.2, 0.3))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$U_{ref}$ time-history')
#plt.xlabel('t/T')
plt.ylabel('$U_{ref}$')

ax2_1.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2_1.grid(b=True, which='major', color='k', linestyle='-')
ax2_1.grid(b=True, which='minor', color='k', linestyle='--')

ax2_2 = fig2.add_subplot(2,1,2)

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', family='Arial')
line_TSR = ax2_2.plot(data2[:,0],data2[:,2],'--',label='$TSR$')
line_TSRmean = ax2_2.plot(mean2[:,0], mean2[:,2],'-',label=r'$\left<TSR\right>_{t}$')
ax2_2.legend(bbox_to_anchor=(0.2, 0.95))
#ax1.legend(handles=[line_thrust,line_power])
#plt.axis([0,40,0,1])
#plt.title(r'$TSR$ time-history')
plt.xlabel('$t$')
plt.ylabel('$TSR$')

ax2_2.get_yaxis().set_minor_locator(mpl.ticker.AutoMinorLocator())
ax2_2.grid(b=True, which='major', color='k', linestyle='-')
ax2_2.grid(b=True, which='minor', color='k', linestyle='--')

plt.savefig('fig2_uref_tsr_'+foldername+'.png')
plt.savefig('fig2_uref_tsr_'+foldername+'.pdf', format='pdf')
#plt.show()
plt.close()


## plot TSR and CT, CP
fig3 = plt.figure()
ax3_1 = fig3.add_subplot(1,1,1)

line_ct = ax3_1.plot(mean2[10000:,2], mean1[10000:,3],'.')
plt.xlabel("$\lambda$")
plt.ylabel("$C_T$")
plt.savefig('fig3_TSR_CT_'+foldername+'.png')
