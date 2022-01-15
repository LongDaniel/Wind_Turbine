import numpy as np
import matplotlib as mpl
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

ti_start = 2600
ti_end = 35000
ti_interval = 200

n_xlist = 3

n_ti = (ti_end - ti_start) / ti_interval + 1

for i in range(1):
  ti = ti_start + i * ti_interval
  fname = "POST_"+"{:010d}".format(ti)+"_UVW_"+"{:04d}".format(3)+".DAT"
  #print(fname)
  data0 = np.genfromtxt(fname)

size = data0.shape
nz = size[0]
nvar = size[1]
print("nz="+str(nz)+", nvar="+str(nvar))

data_all = np.zeros((nz, nvar, n_ti))
time = np.zeros(n_ti)

for i in range(n_ti):
  ti = ti_start + i * ti_interval
  time[i] = ti
  fname = "POST_"+"{:010d}".format(ti)+"_UVW_"+"{:04d}".format(3)+".DAT"
  #print(fname)
  data0 = np.genfromtxt(fname)
  #print(data0)
  data_all[:,:,i] = data0

vel_mean = np.zeros((nz,3))
for i in range(nz):
  for j in range(3):
    vel_mean[i,j] = np.average(data_all[i,3+j,:])

#print("vel_mean:")
#print(vel_mean)


fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)

line_transient = ax1.plot(time, data_all[10,3,:])

plt.savefig('fig_u_timehistory.png')
