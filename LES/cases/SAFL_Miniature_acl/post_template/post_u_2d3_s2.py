import numpy as np
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

import tecplot_io as tec

import post_mplconfig
import os

foldername = os.path.relpath(".","..")

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

nx = 512 
ny = 1
nz = 65
nvar = 10 

#x, y, z, u, v, w = tec.tecplot_reader('umean_2d.dat', 6, [65, 1, 256], 2)
data0 = tec.tecplot_reader('umean_2d.dat', [nz, ny, nx, nvar], 2)
#print(x)
#print(x.reshape(65,256))

x=data0[0,0,:,0].reshape(nx,-1)
z=data0[:,0,0,2].reshape(nz,-1)
u=data0[:,0,:,3].reshape(nz,nx)
v=data0[:,0,:,4].reshape(nz,nx)
w=data0[:,0,:,5].reshape(nz,nx)

xold = x.copy()
xstart0 = x[0]

print(x)

nshift = 0 
idx1 = np.arange(nshift, nx)
idx2 = np.arange(0, nshift)
idx3 = np.concatenate((idx1,idx2))
#print(idx3)

x = x[idx3]
u = u[:, idx3]
v = v[:, idx3]
w = w[:, idx3]
#print(x)
#print(z)
#print(u)

#print(z[27:29])
#u1 = (u[27,31])
#print(u1)
#u_scale =  (2.2/u_infty) / u1


xstart1 = x[0]

xnew = x.copy()
x = xold.copy()
xshift = xstart1 - xstart0

## Mean vs yz
fig_width_pt = 550.0
hratio = 0.25
inches_per_pt = 1.0/72.27
fig_width = fig_width_pt * inches_per_pt
fig_height = fig_width * hratio
#figsize=[fig_width, fig_height]

fig6 = plt.figure(figsize=[fig_width, fig_height])

xtemp = (x-(0.375-xshift)-4.32)/0.15
print(xtemp)
ztemp = z/0.15
X, Z = np.meshgrid(xtemp, ztemp, indexing='ij')
#mean_2d = np.transpose(u) * 2.54390548295 / 2.2
mean_2d = np.transpose(u)

# range of rotor
zrange=np.arange(13,40)
print(z[zrange])
mean_1d = np.zeros((nx))
for i in range(nx):
  mean_1d[i] = np.average(mean_2d[i,zrange]) 

u1 = mean_1d[31]
print(u1)
u_scale =  (2.2/u_infty) / u1

mean_2d = mean_2d * u_scale * u_infty
mean_1d = mean_1d * u_scale

#print(x.shape)
#print(z.shape)
#print(X.shape)
#print(Z.shape)
#print(u.shape)
print(mean_2d.shape)

#ax6_1 = fig6.add_subplot(1,1,1)
#levels = np.arange(0.5, 1.4, 0.1) 
levels = np.arange(1.0, 2.58, 0.04) 
#contour_umean_2d = ax6_1.contourf(X, Z, mean_2d, 10, cmap=plt.cm.YlOrRd, label='U')
plt.contourf(X, Z, mean_2d, levels, cmap=plt.cm.jet, label='U')
plt.colorbar()
#print(plt.figsize)
#plt.figsize=[fig_width, fig_height]
#print(plt.figsize)
plt.xlim([-2.0, 20.0])
plt.ylim([0.0, 2.0])
#ax6_1.legend()
#plt.title('Mean velocity contour')
plt.xlabel('$x/D$')
plt.ylabel('$z/D$')
#plt.axes().set_aspect(1.0, 'datalim')
plt.tight_layout()

#um_hub = np.interp(0.125, z, data_umean)
#print("Um_hub="+str(um_hub))

plt.savefig('fig_6_mean_2d_'+foldername+'.png')
plt.savefig('fig_6_mean_2d_'+foldername+'.pdf', format='pdf')


## Mean vs x
# range of rotor
zrange=np.arange(13,40)
print(z[zrange])
mean_1d = np.zeros((nx))
for i in range(nx):
  mean_1d[i] = np.average(mean_2d[i,zrange]) / 2.2

fig8 = plt.figure()
plt.plot(xtemp, mean_1d)
plt.ylim([0.6, 1.1])
plt.xlabel('$x/D$')
plt.ylabel(r'$U_{rel}=\left<U\right>/U_{hub}$')

i1 = 23 
print(xtemp[i1])
str1 = '${:.3f}'.format(float(mean_1d[i1]))+'\ (x=-2D)$'
plt.annotate(str1, xy=(xtemp[i1], mean_1d[i1]), xytext=(xtemp[i1]-2.0,\
  mean_1d[i1]+0.05), arrowprops=dict(arrowstyle='->'),)

i2 = 41 
print(xtemp[i2])
str1 = '${:.3f}'.format(float(mean_1d[i2]))+'\ (x=0)$'
plt.annotate(str1, xy=(xtemp[i2], mean_1d[i2]), xytext=(xtemp[i2]-6.0,\
  mean_1d[i2]-0.05), arrowprops=dict(arrowstyle='->'),)

i3 = 59 
print(xtemp[i3])
str1 = '${:.3f}'.format(float(mean_1d[i3]))+'\ (x=2D)$'
plt.annotate(str1, xy=(xtemp[i3], mean_1d[i3]), xytext=(xtemp[i3]+2.5,\
  mean_1d[i3]), arrowprops=dict(arrowstyle='->'),)

i4 = 86 
print(xtemp[i4])
str1 = '${:.3f}'.format(float(mean_1d[i4]))+'\ (x=5D)$'
plt.annotate(str1, xy=(xtemp[i4], mean_1d[i4]), xytext=(xtemp[i4]+1.5,\
  mean_1d[i4]-0.05), arrowprops=dict(arrowstyle='->'),)

plt.savefig('fig_8_mean_1d_'+foldername+'.png')
plt.savefig('fig_8_mean_1d_'+foldername+'.pdf', format='pdf')
