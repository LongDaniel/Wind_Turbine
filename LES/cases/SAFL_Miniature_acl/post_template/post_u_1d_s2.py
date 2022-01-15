import numpy as np
from scipy.optimize import curve_fit
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig
import os

def func1(x, a, b):
	return a*x+b

foldername = os.path.relpath(".","..")

## basic parameters
hbar = 0.46
u_infty = 2.54390548295
u_tau = 0.102  # what the others claimed
u_tau1 = 0.102950230486 # what is imported at the inlet
print('u_tau='+str(u_tau))
NU = 1.511e-5
d_rotor = 0.15
um_hub_0 = 2.0251609 / u_infty
um_hub = 2.0251609
um_hub1 = 0.8081174890878188 * u_infty # Used to scale mean velocity in simulation
#um_hub = 0.761218328913 * u_infty 
print("Um_hub="+str(um_hub_0)+"*"+str(u_infty)+"="+str(um_hub))

#u_tau=0.0824319595338
#Um_hub=0.761218328913*2.5439055=1.93646749362

## read experiment result
data_exp1 = np.genfromtxt('../../GetData/fig1_6_exp_-1d.dat', skip_header=4)
data_exp2 = np.genfromtxt('../../GetData/fig1_1_exp_2d.dat', skip_header=4)
data_exp3 = np.genfromtxt('../../GetData/fig1_3_exp_5d.dat', skip_header=4)
data_exp4 = np.genfromtxt('../../GetData/fig2_1_exp_2d.dat', skip_header=4)
data_exp5 = np.genfromtxt('../../GetData/fig2_3_exp_5d.dat', skip_header=4)
data_exp6 = np.genfromtxt('../../GetData/fig3_1_exp_2d.dat', skip_header=4)
data_exp7 = np.genfromtxt('../../GetData/fig3_3_exp_5d.dat', skip_header=4)

## read post-process results 
with np.load('./POST_U_1D1_0004/result_mean.npz') as data:
  z1 = data['arr_0']
  data_umean1 = data['arr_1'] * u_infty 
  data_vmean1 = data['arr_2'] * u_infty
  data_wmean1 = data['arr_3'] * u_infty
  data_uu_mean1 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0008/result_mean.npz') as data:
  z2 = data['arr_0']
  data_umean2 = data['arr_1'] * u_infty 
  data_vmean2 = data['arr_2'] * u_infty
  data_wmean2 = data['arr_3'] * u_infty
  data_uu_mean2 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0011/result_mean.npz') as data:
  z3 = data['arr_0']
  data_umean3 = data['arr_1'] * u_infty 
  data_vmean3 = data['arr_2'] * u_infty
  data_wmean3 = data['arr_3'] * u_infty
  data_uu_mean3 = data['arr_4'] * u_infty**2

_z = z1
_n = len(_z)
i_min = 0
i_max = _n
for i in range(_n):
	if _z[i]<(0.25*d_rotor) and i>i_min :
		i_min = i
	if _z[i]>(2.0*d_rotor) and i<i_max :
		i_max = i
print('i_min='+str(i_min)+', i_max='+str(i_max))

yfit = -data_uu_mean1[i_min:i_max,0,2]
xfit = z1[i_min:i_max]/d_rotor
xfit2 = z1 / d_rotor
popt1, pcov1 = curve_fit(func1, xfit, yfit)
yfit2 = func1(xfit2, *popt1)
u_tau2 = np.sqrt(func1(0, *popt1))
print('u_tau2 = '+str(u_tau2))

## do some adjustment, like non-dimensionalisation
z0 = 0.125
uz0 = np.interp(z0/d_rotor, data_exp1[:,1], data_exp1[:,0])
uz1 = np.interp(z0, z1, data_umean1)
uz2 = np.interp(z0, z2, data_umean2)
uz3 = np.interp(z0, z3, data_umean3)

scale = uz0/(uz1/um_hub)

# time=...

#time=np.arrange(35001,100001)

## fig7: Umean in diameter (no adjustment)

fig7 = plt.figure()

ax1 = fig7.add_subplot(1,2,1)
line1 = ax1.plot(data_umean1/um_hub1, z1/d_rotor, '--')
line2 = ax1.plot(data_exp1[:,0], data_exp1[:,1], '^')
line3 = ax1.plot(data_umean2/um_hub1, z2/d_rotor, '-')
line4 = ax1.plot(data_exp2[:,0], data_exp2[:,1], 'o')

plt.annotate('$x/D=2$', xy=(0.2, 1.75))

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

ax2 = fig7.add_subplot(1,2,2)
line1 = ax2.plot(data_umean1/um_hub1, z1/d_rotor, '--')
line2 = ax2.plot(data_exp1[:,0], data_exp1[:,1], '^')
line3 = ax2.plot(data_umean3/um_hub1, z3/d_rotor, '-')
line4 = ax2.plot(data_exp3[:,0], data_exp3[:,1], 'o')

plt.annotate('$x/D=5$', xy=(0.2, 1.75))

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

plt.savefig('fig_7_umean2_'+foldername+'_original.pdf', format='pdf')
plt.savefig('fig_7_umean2_'+foldername+'_original.png')

## fig8: Umean in diameter (with adjustment)

fig8 = plt.figure()

ax3 = fig8.add_subplot(1,2,1)
line1 = ax3.plot(data_umean1/um_hub1*scale, z1/d_rotor, '--', label='sim., x=-1D')
line2 = ax3.plot(data_exp1[:,0], data_exp1[:,1], '^', label='exp., x=-1D')
line3 = ax3.plot(data_umean2/um_hub1*scale, z2/d_rotor, '-', label='sim., x=2D')
line4 = ax3.plot(data_exp2[:,0], data_exp2[:,1], 'o', label='exp., x=2D')
ax3.legend(loc='upper left')
plt.annotate('$x/D=2$', xy=(1.1, 0.2))

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

ax4 = fig8.add_subplot(1,2,2)
line1 = ax4.plot(data_umean1/um_hub1*scale, z1/d_rotor, '--', label='sim., x=-1D')
line2 = ax4.plot(data_exp1[:,0], data_exp1[:,1], '^', label='exp., x=-1D')
line3 = ax4.plot(data_umean3/um_hub1*scale, z3/d_rotor, '-', label='sim., x=5D')
line4 = ax4.plot(data_exp3[:,0], data_exp3[:,1], 'o', label='exp., x=5D')
ax4.legend(loc='upper left')
plt.annotate('$x/D=5$', xy=(1.1, 0.2))

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

plt.savefig('fig_7_umean2_'+foldername+'.pdf', format='pdf')
plt.savefig('fig_7_umean2_'+foldername+'.png')


## fig9: u_rms
fig9 = plt.figure()

ax5_1 = fig9.add_subplot(1,2,1)
line1 = ax5_1.plot(np.sqrt(data_uu_mean2[:,0,0])/u_tau2, z2/d_rotor,'-')
line2 = ax5_1.plot(data_exp4[:,0], data_exp4[:,1],'o')
#line3 = ax5_1.plot(np.sqrt(data_uu_mean1[:,0,0])/u_tau, z2/d_rotor,'--')
plt.ylim([0,2.0])
plt.xlabel('$u_{rms}^\prime/u_*$')
plt.ylabel('$z/D$')

ax5_2 = fig9.add_subplot(1,2,2)
line4 = ax5_2.plot(np.sqrt(data_uu_mean3[:,0,0])/u_tau2, z3/d_rotor,'-')
line5 = ax5_2.plot(data_exp5[:,0], data_exp5[:,1],'o')
#line6 = ax5_2.plot(np.sqrt(data_uu_mean1[:,0,0])/u_tau, z3/d_rotor,'--')
plt.ylim([0,2.0])
plt.xlabel('$u_{rms}^\prime/u_*$')
plt.ylabel('$z/D$')
plt.savefig('fig_9_urms_'+foldername+'.png')

## savedata
np.savez('result_umean_z',\
  data_umean1/um_hub1*scale,z1/d_rotor,data_exp1[:,0],data_exp1[:,1],\
  data_umean2/um_hub1*scale,z2/d_rotor,data_exp2[:,0],data_exp2[:,1],\
  data_umean3/um_hub1*scale,z3/d_rotor,data_exp3[:,0],data_exp3[:,1])
np.savez('result_urms_z',\
  np.sqrt(data_uu_mean2[:,0,0])/u_tau2, z2/d_rotor,data_exp4[:,0], data_exp4[:,1],\
  np.sqrt(data_uu_mean3[:,0,0])/u_tau2, z3/d_rotor,data_exp5[:,0], data_exp5[:,1])

## fig11: u_rms
fig11 = plt.figure()

ax7_1 = fig11.add_subplot(1,2,1)
line1 = ax7_1.plot(np.sqrt(data_uu_mean2[:,0,0])/um_hub1, z2/d_rotor,'-', label='sim., x=2D')
line2 = ax7_1.plot(data_exp4[:,0]*u_tau2/um_hub1, data_exp4[:,1],'o', label='exp., x=2D')
line3 = ax7_1.plot(np.sqrt(data_uu_mean1[:,0,0])/um_hub1, z2/d_rotor,'--', label='sim., x=-1D')
ax7_1.legend(loc='lower left')
plt.ylim([0,2.0])
plt.xlabel('$u_{rms}^\prime/U_{hub}$')
plt.ylabel('$z/D$')

ax7_2 = fig11.add_subplot(1,2,2)
line4 = ax7_2.plot(np.sqrt(data_uu_mean3[:,0,0])/um_hub1, z3/d_rotor,'-', label='sim., x=5D')
line5 = ax7_2.plot(data_exp5[:,0]*u_tau2/um_hub1, data_exp5[:,1],'o', label='exp., x=5D')
line6 = ax7_2.plot(np.sqrt(data_uu_mean1[:,0,0])/um_hub1, z3/d_rotor,'--', label='sim, x=-1D')
ax7_2.legend(loc='lower left')
plt.ylim([0,2.0])
plt.xlabel('$u_{rms}^\prime/U_{hub}$')
plt.ylabel('$z/D$')
plt.savefig('fig_11_urms_'+foldername+'.png')

## fig10: u'w'
fig10 = plt.figure()
ax6_1 = fig10.add_subplot(1,2,1)
line1 = ax6_1.plot(-data_uu_mean2[:,0,2]/u_tau2/u_tau2, z2/d_rotor,'-')
line2 = ax6_1.plot(data_exp6[:,0], data_exp6[:,1],'o')
line3 = ax6_1.plot(-data_uu_mean1[:,0,2]/u_tau2/u_tau2, z1/d_rotor,'--')
line4 = ax6_1.plot(yfit2/u_tau2/u_tau2, xfit2,':')
plt.ylim([0,2.0])
plt.xlabel(r'$-\overline{u^\prime w^\prime}/u_*^2$')
plt.ylabel('$z/D$')

ax6_2 = fig10.add_subplot(1,2,2)
plt.ylim([0,2.0])
line3 = ax6_2.plot(-data_uu_mean3[:,0,2]/u_tau2/u_tau2, z3/d_rotor,'-')
line4 = ax6_2.plot(data_exp7[:,0], data_exp7[:,1],'o')
plt.xlabel(r'$-\overline{u^\prime w^\prime}/u_*^2$')
plt.ylabel('$z/D$')
plt.savefig('fig_10_uw_'+foldername+'.png')
