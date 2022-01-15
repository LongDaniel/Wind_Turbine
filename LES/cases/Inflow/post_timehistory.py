import numpy as np
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

## basic parameters
hbar = 0.46
u_infty = 2.5439055
u_tau = np.sqrt(0.002*0.96) * u_infty
print('u_tau='+str(u_tau))
RE = 168359.066
NU = 1.511e-5

## read post-process results 
with np.load('result_timehistory.npz') as data:
#  np.savez('result_timehistory', z, data_u_output, iz, iz_array)
  z = data['arr_0']
  data_u_output = data['arr_1'] * u_infty 
  iz = data['arr_2'] 
  iz_array = data['arr_3'] 

time=np.arange(35001,100001) * 0.001

## fig1: Instantaneous vs z,t
fig1 = plt.figure()

#ax1_1 = fig1.add_subplot(3,1,1)
#line_u_z1 = ax1_1.plot(time, data_u[:,iz-4])

#ax1_2 = fig1.add_subplot(3,1,2)
#line_u_z2 = ax1_2.plot(time, data_u[:,iz])

#ax1_3 = fig1.add_subplot(3,1,3)
#line_u_z3 = ax1_3.plot(time, data_u[:,iz+5])

ax1_1 = fig1.add_subplot(1,1,1)
line_u_z1 = ax1_1.plot(time, data_u_output[:, 0]/u_infty,label='$z/H='+'{:0.3f}'.format(z[iz_array[0]]/hbar)+'$')
line_u_z2 = ax1_1.plot(time, data_u_output[:, 1]/u_infty,label='$z/H='+'{:0.3f}'.format(z[iz_array[1]]/hbar)+'$')
line_u_z3 = ax1_1.plot(time, data_u_output[:, 2]/u_infty,label='$z/H='+'{:0.3f}'.format(z[iz_array[2]]/hbar)+'$')
ax1_1.legend()
#plt.title('Streamwise velocity at various heights')
plt.xlabel('$t$')
plt.ylabel('$u/U_\infty$')

plt.savefig('fig_1_instantaneous_various_z.png')
plt.savefig('fig_1_instantaneous_various_z.pdf', format='pdf')


## fig4: quadrant
#fig4 = plt.figure()
#ax4_1 = fig4.add_subplot(1,1,1)
#line_u_w = ax4_1.plot(data_velturb_output[:,:,0], data_velturb_output[:,:,2],'b.')
#plt.xlabel("u'")
#plt.ylabel("w'")
#plt.savefig('fig_4_quadrant.png')

