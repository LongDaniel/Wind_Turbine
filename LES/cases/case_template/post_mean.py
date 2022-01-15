import sys
import h5py
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def func2(x, a, b):
  return a*x+b

def find_index(_z, _limits):
  _n = len(_z)
  _i_min = 0
  _i_max = _n - 1
  _limits2 = np.zeros(2)
  if isinstance(_limits, float):
    _limits2[0:2] = _limits
  else:
    _limits2 = _limits
          
  for i in range(_n):
    if _z[i]<_limits2[0] and i>_i_min :
      _i_min = i
    if _z[i]>_limits2[1] and i<_i_max :
      _i_max = i
  #print('zlimits='+str(_limits))
  #print('i_min='+str(_i_min)+', i_max='+str(_i_max))
  
  if isinstance(_limits, float):
    return _i_min
  else:
    return _i_min, _i_max

def diff_central(x, y):
  x0 = x[:-2]
  x1 = x[1:-1]
  x2 = x[2:]
  y0 = y[:-2]
  y1 = y[1:-1]
  y2 = y[2:]
  f = (x2 - x1)/(x2 - x0)
  f1 = (1-f)*(y2 - y1)/(x2 - x1) + f*(y1 - y0)/(x1 - x0)
  f2 = x.copy()
  f2[1:-1] = f1
  f2[0] = f1[0]
  f2[-1] = f1[-1]
  return f2
  
# for two-column layout, the width is 3.487 inches
fig_width_pt = 800.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_myratio = 0.85
fig_width = fig_width_pt*inches_per_pt  # width in inches
#fig_height = fig_width      # height in inches
#fig_height = fig_width*golden_mean      # height in inches
#fig_height = fig_width/golden_mean      # height in inches
fig_height = fig_width * fig_myratio
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps',
          'font.size': 18,
          'axes.labelsize': 22,
          'text.fontsize': 22,
          'legend.fontsize': 14,
          'xtick.labelsize': 16,
          'ytick.labelsize': 16,
          'text.usetex': True,
          'figure.figsize': fig_size
#          'legend.unicode_minus': False
}
#          'lines.markerfacecolor': 'none',
#          'scatter.markerfacecolor': 'none'
mpl.rcParams.update(params)

#NPX = 256
#NPY = 64
#NPZ = 65

kappa = 0.4
nu = 1.511e-5

PEX = 1.45444104333
PEY = 8.72664625997
hbar = 0.46

uinfty = 2.54390548295

mode = 1  

if mode==1:
  g = h5py.File('grid.h5', 'r')
  zw = g["zw"]

  tis = 40100 
  #tis = 100000
  tie = 50000
  tii = 100
  nti = int((tie - tis) / tii + 1)

  for it in range(nti):
    ti = tis + tii * it
    fname = 'DAT_{:010d}.h5'.format(ti)
    print("Reading file "+ fname)
    f = h5py.File(fname, "r")
    
    #print("Keys: %s" % f.keys())
    ## Old version Keys: [u'dz', u'dzw', u'eta', u'eta0', u'hh', u'pp', u'u', u'v', u'w', u'z', u'zw', u'zz']
    ## New version Keys: [u'eta', u'hh', u'pp', u'u', u'v', u'w', u'z']

    zz = np.array(f["z"][:,0,0]).copy()
    ##zw = f["zw"]
    u = f["u"]
    v = f["v"]
    w2 = f["w"]
    w = np.array(w2).copy()

    print(u.shape)
    ##(65, 64, 256)

    #print(np.vstack((zz, zw)))

    NPX = u.shape[2]
    NPY = u.shape[1]
    NPZ = u.shape[0]
    
    if it==0:
      u_m_all = np.zeros((NPZ,3))
      uu_m_all = np.zeros((NPZ,6))

    w[0, :, :] = w2[0, :, :]
    for k in range(1,NPZ):
      w[k, :, :] = 0.5*(w2[k-1, :, :] + w2[k, :, :])

    u_m = np.zeros((NPZ,3))
    uu_m1 = np.zeros((NPZ,6))
    uu_m2 = np.zeros((NPZ,6))

    for k in range(NPZ):
      u_m[k,0] = np.average(u[k,:,:])
      u_m[k,1] = np.average(v[k,:,:])
      u_m[k,2] = np.average(w[k,:,:])
      uu_m1[k,0] = np.average(u[k,:,:]**2)
      uu_m1[k,1] = np.average(v[k,:,:]**2)
      uu_m1[k,2] = np.average(w[k,:,:]**2)
      uu_m1[k,3] = np.average(u[k,:,:]*v[k,:,:])
      uu_m1[k,4] = np.average(u[k,:,:]*w[k,:,:])
      uu_m1[k,5] = np.average(v[k,:,:]*w[k,:,:])
      uu_m2[k, 0:3] = uu_m1[k, 0:3] - u_m[k,0:3]**2
      uu_m2[k, 3] = uu_m1[k, 3] - u_m[k, 0] * u_m[k, 1]
      uu_m2[k, 4] = uu_m1[k, 4] - u_m[k, 0] * u_m[k, 2]
      uu_m2[k, 5] = uu_m1[k, 5] - u_m[k, 1] * u_m[k, 2]
      
    u_m_all = u_m_all + u_m
    uu_m_all = uu_m_all + uu_m2
    f.close()

  u_m_all = u_m_all / nti
  uu_m_all = uu_m_all / nti

  fig1 = plt.figure()
  ax = fig1.add_subplot(1,2,1)
  labels=["U", "V", "W"]
  for i in range(3):
    ax.plot(u_m_all[:,i], zz, label=labels[i])
  plt.legend()

  ax = fig1.add_subplot(1,2,2)
  labels=["uu", "vv", "ww", "uv", "uw", "vw"]
  for i in range(6):
    ax.plot(uu_m_all[:,i], zz, label=labels[i])
  plt.legend()

  fig1.savefig("mean.png")
  plt.close()

  g.close()

  np.savez("umean", u_m_all=u_m_all, uu_m_all=uu_m_all, zz=zz, NPX=NPX, NPY=NPY, NPZ=NPZ)
elif mode==2:
  npzfiles = np.load("umean.npz")
  u_m_all = npzfiles["u_m_all"]
  uu_m_all = npzfiles["uu_m_all"]
  zz = npzfiles["zz"] 

  u_m_all = u_m_all * uinfty
  uu_m_all = uu_m_all * uinfty**2
  
  # fitting <u'w'>
  imin1, imax1 = find_index(zz, np.array((0.125,0.35)))
  yfit1 = uu_m_all[imin1:imax1, 4]
  xfit1 = zz[imin1:imax1]
  popt1, pcov1 = curve_fit(func2, xfit1, yfit1)
  xfitout1 = zz
  yfitout1 = func2(xfitout1, *popt1)
  ufric = np.sqrt(-yfitout1[0])
  print("ufric="+str(ufric)+", ufric/uinfty="+str(ufric/uinfty))  
  
  zp = (zz) * (ufric) / nu
  dudz = diff_central(zz, u_m_all[:, 0])
  phi = zz / ufric * dudz 
  
  print("i\tzz\tz+\tum\tu+")
  for i in range(len(zz)):
    print("{:d}\t{:0.4f}\t{:0.4f}\t{:0.4f}\t{:0.4f}".format(i, zz[i], zp[i], u_m_all[i,0], u_m_all[i,0]/ufric))
    #print(str(i)+"\t"+str(zz[i])+"\t"+str(zp[i])+"\t"+str(u_m_all[i,0])) 

  fig1 = plt.figure()
  ax = fig1.add_subplot(2,2,1)
  labels=["U", "V", "W"]
  for i in range(3):
    ax.plot(u_m_all[:,i], zz, label=labels[i])
  plt.legend()
  
  ax = fig1.add_subplot(2,2,2)
  labels=["uu", "vv", "ww", "uv", "uw", "vw"]
  for i in range(6):
    ax.plot(uu_m_all[:,i]/(ufric**2), zz, label=labels[i])
  plt.legend()
  
  ax = fig1.add_subplot(2,2,3)
  ax.plot(zp, u_m_all[:,0]/ufric, '-+')
  plt.xscale('log')
  
  ax = fig1.add_subplot(2,2,4)
  ax.plot(phi * kappa, zp, '-+')
  plt.ylim([0.0, 200])

  fig1.savefig("mean_load.png")
  plt.close()
  
  wallu = np.vstack((zz, u_m_all[:,0], zp, u_m_all[:,0]/ufric)).transpose()
  np.savetxt("wallu.dat", wallu)
  
  Tt = -uu_m_all[:,4]
  Tv = nu * dudz
  T0 = ufric**2
  
  Ttp = Tt / T0
  Tvp = Tv / T0
  phi2 = kappa * zp *(1-Ttp-zz/hbar)
  
  fig2 = plt.figure()
  ax = fig2.add_subplot(1,2,1)
  ax.plot(Tt/T0, zp, '-o')
  ax.plot(Tv/T0, zp, '-.')
  ax.plot(Tt/T0+Tv/T0, zp, '--')
  plt.ylim([0.0, 200])
  
  
  ax = fig2.add_subplot(1,2,2)
  ax.plot(phi2, zp, '-+')
  ax.plot(phi, zp, '-.')
  plt.ylim([0.0, 200])
  
  fig2.savefig("phi_load.png")
  plt.close()
    
  
