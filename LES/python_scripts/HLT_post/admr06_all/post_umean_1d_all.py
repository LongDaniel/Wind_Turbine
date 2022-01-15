import numpy as np
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

## basic parameters
hbar = 0.46
u_infty = 2.5439055
u_tau = np.sqrt(0.00105) * u_infty
print('u_tau='+str(u_tau))
RE = 168359.066
NU = 1.511e-5
d_rotor = 0.15
um_hub_0 = 0.761218328913 
um_hub = 2.2
#um_hub = 0.761218328913 * u_infty 
print("Um_hub="+str(um_hub_0)+"*"+str(u_infty)+"="+str(um_hub))

#u_tau=0.0824319595338
#Um_hub=0.761218328913*2.5439055=1.93646749362

## read experiment result
data_exp1 = np.genfromtxt('./GetData/fig1_6_exp_-1d.dat', skip_header=4)
data_exp2 = np.genfromtxt('./GetData/fig1_1_exp_2d.dat', skip_header=4)
data_exp3 = np.genfromtxt('./GetData/fig1_3_exp_5d.dat', skip_header=4)

## read post-process results 
with np.load('./POST_U_1D1_0001/result_mean.npz') as data:
  z1 = data['arr_0']
  data_umean1 = data['arr_1'] * u_infty 
  data_vmean1 = data['arr_2'] * u_infty
  data_wmean1 = data['arr_3'] * u_infty
  data_uu_mean1 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0005/result_mean.npz') as data:
  z2 = data['arr_0']
  data_umean2 = data['arr_1'] * u_infty 
  data_vmean2 = data['arr_2'] * u_infty
  data_wmean2 = data['arr_3'] * u_infty
  data_uu_mean2 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0008/result_mean.npz') as data:
  z3 = data['arr_0']
  data_umean3 = data['arr_1'] * u_infty 
  data_vmean3 = data['arr_2'] * u_infty
  data_wmean3 = data['arr_3'] * u_infty
  data_uu_mean3 = data['arr_4'] * u_infty**2
# time=...

#time=np.arrange(35001,100001)

## fig7: Umean in diameter

fig7 = plt.figure()

ax1 = fig7.add_subplot(1,2,1)
line1 = ax1.plot(data_umean1/um_hub, z1/d_rotor, '--')
line2 = ax1.plot(data_exp1[:,0], data_exp1[:,1], '^')
line3 = ax1.plot(data_umean2/um_hub, z2/d_rotor, '-')
line4 = ax1.plot(data_exp2[:,0], data_exp2[:,1], 'o')

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

ax2 = fig7.add_subplot(1,2,2)
line1 = ax2.plot(data_umean1/um_hub, z1/d_rotor, '--')
line2 = ax2.plot(data_exp1[:,0], data_exp1[:,1], '^')
line3 = ax2.plot(data_umean3/um_hub, z3/d_rotor, '-')
line4 = ax2.plot(data_exp3[:,0], data_exp3[:,1], 'o')

plt.xlim([0,1.5])
plt.ylim([0,2.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

plt.savefig('fig_7_umean2.pdf', format='pdf')
plt.savefig('fig_7_umean2.png')


 
