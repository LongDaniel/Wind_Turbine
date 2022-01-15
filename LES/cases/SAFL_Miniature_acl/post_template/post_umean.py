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
um_hub = 0.761218328913 * u_infty 
print("Um_hub="+str(um_hub_0)+"*"+str(u_infty)+"="+str(um_hub))

#u_tau=0.0824319595338
#Um_hub=0.761218328913*2.5439055=1.93646749362

## read experiment result
data_exp = np.genfromtxt('../GetData/fig1_6_exp_-1d.dat', skip_header=4)

## read post-process results 
with np.load('result_mean.npz') as data:
  z = data['arr_0']
  data_umean = data['arr_1'] * u_infty 
  data_vmean = data['arr_2'] * u_infty
  data_wmean = data['arr_3'] * u_infty
  data_uu_mean = data['arr_4'] * u_infty**2
# time=...

#time=np.arrange(35001,100001)

## fig3: 2nd order moments vs. z
fig3 = plt.figure()

ax3_1 = fig3.add_subplot(1,1,1)
#plt.rc('text', usetex=True)
#plt.rc('font', family='Arial')

line_uu_z = ax3_1.plot(data_uu_mean[:,0,0]/(u_tau**2), z/hbar,'-', label="$\overline{u^\prime u^\prime}/u_*^2$")
line_vv_z = ax3_1.plot(data_uu_mean[:,1,1]/(u_tau**2), z/hbar,':', label="$\overline{v^\prime v^\prime}/u_*^2$")
line_ww_z = ax3_1.plot(data_uu_mean[:,2,2]/(u_tau**2), z/hbar,'-.', label="$\overline{w^\prime w^\prime}/u_*^2$")
line_uw_z = ax3_1.plot(-data_uu_mean[:,0,2]/(u_tau**2), z/hbar,'--', label="$-\overline{u^\prime w^\prime}/u_*^2$")
#line_vv_z = ax3_1.plot(data_uu_mean[:,1,1], z, label=r"$\overline{v^\prime v^\prime}$")
#line_ww_z = ax3_1.plot(data_uu_mean[:,2,2], z, label=r"$\overline{w^\prime w^\prime}$")
#line_uw_z = ax3_1.plot(-data_uu_mean[:,0,2], z, label=r"$\overline{-u^\prime w^\prime}$")

#line_uu_z = ax3_1.plot(data_uu_mean[:,0,0], z, label=r"$u'u'$")
#line_vv_z = ax3_1.plot(data_uu_mean[:,1,1], z, label=r"$v'v'$")
#line_ww_z = ax3_1.plot(data_uu_mean[:,2,2], z, label=r"$w'w'$")
#line_uw_z = ax3_1.plot(data_uu_mean[:,0,2], z, label=r"$u'w'$")

ax3_1.legend()
#plt.title('Second order velocity moments vs. height')
plt.xlabel('$\overline{u_i^\prime u_j^\prime}/u_*^2$')
plt.ylabel('$z/H$')

plt.savefig('fig_3_2ndmoment_z.png')
plt.savefig('fig_3_2ndmoment_z.pdf', format='pdf')


## fig6: Umean in wall coordinate

fig6 = plt.figure()
ax1 = fig6.add_subplot(1,2,1)
line1 = ax1.plot(data_umean, z, 'x-')
#ax1.set_xscale('log')
plt.xlabel('$U(z)$')
plt.ylabel('$z$')

z_u = np.column_stack((z,data_umean))
#print(z_u)

uplus = data_umean / u_tau
yplus = z * u_tau / NU

yplus_uplus = np.column_stack((yplus,uplus))
#print(yplus_uplus)

ax2 = fig6.add_subplot(1,2,2)
line2 = ax2.plot(yplus, uplus)
ax2.set_xscale('log')

#ax1.set_yscale('log')
##plt.show()
#
##print(a['freq'])
##print(a['sp_abs'])
#
#
con1 = (yplus>=200.0)
con2 = (yplus<1000.0)
con3 = con1 & con2
#
con4 = (yplus>2.0)
con5 = (yplus<3000.0)
con6 = con4 & con5
##print(con1)
##print(con2)
##print(con3)
#
b1 = np.where(con3)[0]
b2 = np.where(con6)[0]
##print(b1)
#
yplus2 = yplus[b1]
yplus3 = yplus[b2]
yplus2_log = np.log(yplus2)
uplus2 = uplus[b1]
#
##print(freq2_log)
#
## https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
coeff1 = np.polyfit(yplus2_log, uplus2, 1)
print(coeff1)
#
fun_lin = np.poly1d(coeff1)
fun_lin_expr = '$u^+ = '+'\\dfrac{1}{'+'{:.4f}'.format(1.0/coeff1[0])+'} \ln (y^+) +'+'{:.4f}$'.format(coeff1[1])
expr_sim = '$\\dfrac{1}{\kappa}='+'{:.2f}$'.format(coeff1[0])
#
uplus3 = fun_lin(np.log(yplus3))


line2 = ax2.plot(yplus3, uplus3, '--')
#
#
ax2.annotate(fun_lin_expr, xy=(30, 8))
  
plt.xlabel('$y^+$')
plt.ylabel('$u^+$')

#plt.xlim(0.1,100)
#plt.ylim(-0.01,1000)
plt.savefig('fig_6_wall.pdf', format='pdf')
plt.savefig('fig_6_wall.png')

## fig7: Umean in diameter

fig7 = plt.figure()
ax1 = fig7.add_subplot(1,1,1)
line1 = ax1.plot(data_umean/um_hub, z/d_rotor, 'x-')
line2 = ax1.plot(data_exp[:,0], data_exp[:,1], '^')

plt.ylim([0,2.0])

#ax1.set_xscale('log')
plt.xlabel('$U(z)/U_{hub}$')
plt.ylabel('$z/D$')

plt.savefig('fig_7_umean2.pdf', format='pdf')
plt.savefig('fig_7_umean2.png')


 
