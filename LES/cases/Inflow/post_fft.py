import numpy as np
import matplotlib as mpl

# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import post_mplconfig

with np.load('result_fft.npz') as data:
  freq = data['arr_0']
  sp = data['arr_1']**2

fig1 = plt.figure()
ax1 = fig1.add_subplot(1,1,1)
line1 = ax1.plot(freq, sp)
ax1.set_xscale('log')
ax1.set_yscale('log')
#plt.show()

#print(a['freq'])
#print(a['sp_abs'])


con1 = (freq>=5.0)
con2 = (freq<15.0)
con3 = con1 & con2

con4 = (freq>2.0)
con5 = (freq<60.0)
con6 = con4 & con5
#print(con1)
#print(con2)
#print(con3)

b1 = np.where(con3)[0]
b2 = np.where(con6)[0]
#print(b1)

freq2 = freq[b1]
freq3 = freq[b2]
freq2_log = np.log(freq2)
sp2_log = np.log(sp[b1])

#print(freq2_log)

# https://docs.scipy.org/doc/numpy/reference/generated/numpy.polyfit.html
coeff1 = np.polyfit(freq2_log, sp2_log, 1)
print(coeff1)

fun_lin = np.poly1d(coeff1)
fun_lin_expr = 'log(sp) = '+'{:.4f}'.format(coeff1[0])+' x +'+'{:.4f}'.format(coeff1[1])
expr_sim = '$k='+'{:.2f}'.format(coeff1[0])+'$'

sp3_log = fun_lin(np.log(freq3))
sp3 = np.exp(sp3_log)

line2 = ax1.plot(freq3, sp3, '--')


ax1.annotate(expr_sim, xy=(20, 8))
  
plt.xlabel('$f$')
plt.ylabel('$S(f)$')

plt.xlim(0.1,100)
#plt.ylim(-0.01,1000)
plt.savefig('fig_5_fft.pdf', format='pdf')
plt.savefig('fig_5_fft.png')
