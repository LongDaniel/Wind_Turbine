import numpy as np
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

import post_mplconfig
import os

foldername = os.path.relpath(".","..")

## basic parameters
hbar = 0.46
hhub = 0.125
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
data_umean_exp1 = np.genfromtxt('../../GetData/fig1_6_exp_-1d.dat', skip_header=4)
data_uumean_exp1 = np.genfromtxt('../../GetData/TI_PorteAgel_exp_-1d.dat', skip_header=4)

## read post-process results 
with np.load('./POST_U_1D1_0001/result_mean.npz') as data:
  z1 = data['arr_0']
  data_umean1 = data['arr_1'] * u_infty 
  data_vmean1 = data['arr_2'] * u_infty
  data_wmean1 = data['arr_3'] * u_infty
  data_uu_mean1 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0002/result_mean.npz') as data:
  z2 = data['arr_0']
  data_umean2 = data['arr_1'] * u_infty 
  data_vmean2 = data['arr_2'] * u_infty
  data_wmean2 = data['arr_3'] * u_infty
  data_uu_mean2 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0003/result_mean.npz') as data:
  z3 = data['arr_0']
  data_umean3 = data['arr_1'] * u_infty 
  data_vmean3 = data['arr_2'] * u_infty
  data_wmean3 = data['arr_3'] * u_infty
  data_uu_mean3 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0004/result_mean.npz') as data:
  z4 = data['arr_0']
  data_umean4 = data['arr_1'] * u_infty 
  data_vmean4 = data['arr_2'] * u_infty
  data_wmean4 = data['arr_3'] * u_infty
  data_uu_mean4 = data['arr_4'] * u_infty**2

with np.load('./POST_U_1D1_0005/result_mean.npz') as data:
  z5 = data['arr_0']
  data_umean5 = data['arr_1'] * u_infty 
  data_vmean5 = data['arr_2'] * u_infty
  data_wmean5 = data['arr_3'] * u_infty
  data_uu_mean5 = data['arr_4'] * u_infty**2

## do some adjustment, like non-dimensionalisation
#z0 = 2.0*d_rotor
#uz0 = np.interp(z0/d_rotor, data_exp1[:,1], data_exp1[:,0])
#uz1 = np.interp(z0, z1, data_umean1)
#uz2 = np.interp(z0, z2, data_umean2)
#uz3 = np.interp(z0, z3, data_umean3)

#scale = uz0/(uz1/um_hub)

# time=...

#time=np.arrange(35001,100001)

## fig7: Umean in diameter (no adjustment)

fig7 = plt.figure()

ax1 = fig7.add_subplot(1,2,1)
line1_1 = ax1.plot(data_umean_exp1[:,0], data_umean_exp1[:,1]*d_rotor/hbar, '^', label='exp.')
line1_2 = ax1.plot(data_umean1/um_hub1, z1/hbar, '-', label='sim., x=-5.0D')
line1_3 = ax1.plot(data_umean2/um_hub1, z2/hbar, ':', label='sim., x=-4.0D')
line1_4 = ax1.plot(data_umean3/um_hub1, z3/hbar, '--', label='sim., x=-3.0D')
line1_5 = ax1.plot(data_umean4/um_hub1, z4/hbar, '.', label='sim., x=-2.0D')
line1_6 = ax1.plot(data_umean5/um_hub1, z5/hbar, '-.', label='sim., x=-1.0D')
ax1.legend(loc='upper left')
#plt.annotate('$x/D=2$', xy=(0.2, 1.75))

#plt.xlim([0,1.5])
plt.ylim([0,1.0])
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/\delta$')


ax2 = fig7.add_subplot(1,2,2)
line2_1 = ax2.plot(data_uumean_exp1[:,0], data_uumean_exp1[:,1]*d_rotor/hbar, '^', label='exp.')

# fig_sw: 1, normal; 2, scaled
fig_sw = 1

def func1(x, a, b):
  return a*x + b

for i in range(z1.shape[0]):
  print(str(i)+' '+str(z1[i]))

bound = np.array((0.05, 0.3))
data_temp = data_uumean_exp1[:,1]*d_rotor
#print(data_temp>bound[0])
index = np.transpose(np.nonzero((data_temp>bound[0]) & (data_temp<bound[1])))
index = index.reshape(len(index))
#print(index)
#for i in range(index.shape[0]):
#  print(str(i)+' '+str(index[i])+' '+str(data_temp[i]))
ydata1 = data_uumean_exp1[index,0]
xdata1 = data_uumean_exp1[index,1]*d_rotor/hbar
for i in range(len(xdata1)):
  print(str(i)+' '+str(xdata1[i])+' '+str(ydata1[i]))
popt1, pcov1 = curve_fit(func1, xdata1, ydata1)
yfit1 = func1(xdata1, *popt1)

data_temp = z5
#print(data_temp>bound[0])
index = np.transpose(np.nonzero((data_temp>bound[0]) & (data_temp<bound[1])))
index = index.reshape(len(index))
ydata2 = np.sqrt(data_uu_mean5[index, 0,0])/um_hub1
xdata2 = z5[index]/hbar
popt2, pcov2 = curve_fit(func1, xdata2, ydata2)
yfit2 = func1(xdata2, *popt2)
print(popt1)
print(popt2)

temp1 = func1(hhub, *popt1)
temp2 = func1(hhub, *popt2)
pscale = temp1/temp2
print('pscale='+str(pscale))
 
if fig_sw==1:
  line2_2 = ax2.plot(np.sqrt(data_uu_mean1[:,0,0])/um_hub1, z1/hbar, '-', label='sim., x=-5.0D')
  line2_3 = ax2.plot(np.sqrt(data_uu_mean2[:,0,0])/um_hub1, z2/hbar, ':', label='sim., x=-4.0D')
  line2_4 = ax2.plot(np.sqrt(data_uu_mean3[:,0,0])/um_hub1, z3/hbar, '--', label='sim., x=-3.0D')
  line2_5 = ax2.plot(np.sqrt(data_uu_mean4[:,0,0])/um_hub1, z4/hbar, '.', label='sim., x=-2.0D')
  line2_6 = ax2.plot(np.sqrt(data_uu_mean5[:,0,0])/um_hub1, z5/hbar, '-.', label='sim., x=-1.0D')
else:
  line2_6 = ax2.plot(np.sqrt(data_uu_mean5[:,0,0])/um_hub1*pscale, z5/hbar, 'o', label='sim., x=-1.0D')
  line_cf1 = ax2.plot(yfit1, xdata1, '-', label='fit of exp.')
  line_cf2 = ax2.plot(yfit2*pscale, xdata2, '-', label='fit of sim.')
ax2.legend()
#plt.annotate('$x/D=2$', xy=(0.2, 1.75))

#plt.xlim([0,1.5])
plt.ylim([0,1.0])
plt.xlabel(r'$\sqrt{\overline{u^\prime u^\prime}}/U_{hub}$')
#plt.ylabel('$z/\delta$')

plt.savefig('fig_12_inflow_'+foldername+'.pdf', format='pdf')
plt.savefig('fig_12_inflow_'+foldername+'.png')

