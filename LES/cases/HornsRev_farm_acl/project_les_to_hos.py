import shutil
import h5py
import numpy as np
import matplotlib as mpl

## Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

from scipy import fftpack

fname_in_les = "restart.h5"
fname_in_hos = "restart_hos.h5"
#shutil.move(fname_out1, fname_out2)

f_in_les = h5py.File(fname_in_les, "r")
f_in_hos = h5py.File(fname_in_hos, "r")

print("Keys in f_in_les: %s" % f_in_les.keys())
print("Keys in f_in_hos: %s" % f_in_hos.keys())

eta_in_les = f_in_les["eta"]
hh_in_les = f_in_les["hh"]
pp_in_les = f_in_les["pp"]
u_in_les = f_in_les["u"]
v_in_les = f_in_les["v"]
w_in_les = f_in_les["w"]

eta_in_hos = f_in_hos["eta_hos"]
pa_in_hos = f_in_hos["pa_hos"]
vps_in_hos = f_in_hos["vps_hos"]


print(eta_in_les.shape)
print(eta_in_hos.shape)
print(pp_in_les.shape)
size_les2d = np.array((384, 384))
size_hos2d = np.array((768, 768))
size_les3d = np.array((129, 384, 384))

## mode:
## 1. plot only
## 2. do the correction for restart file.

mode = 2 

## For all modes, firstly we plot it.
fig1 = plt.figure()
ax = fig1.add_subplot(3,1,1)
ax.contourf(eta_in_les)
ax = fig1.add_subplot(3,1,2)
ax.contourf(hh_in_les)
ax = fig1.add_subplot(3,1,3)
ax.contourf(eta_in_hos)

plt.savefig("projection_fig1.png")
plt.close()

if mode==2:
	## we found eta in LES.h5 have a wrong sign,  so correct it.
	fname_out_les = "restart_out.h5"
	fname_out_hos = "restart_hos_out.h5"
	f_out_les = h5py.File(fname_out_les, "w")	
	
	dset = f_out_les.create_dataset("eta", size_les2d, dtype='float64')
	dset[:,:] = -1.0 * np.array(hh_in_les)
	
	dset = f_out_les.create_dataset("hh", size_les2d, dtype='float64')
	dset[:,:] = np.array(hh_in_les)
	
	dset = f_out_les.create_dataset("pp", size_les3d, dtype='float64')
	dset[:,:,:] = np.array(pp_in_les)
	
	dset = f_out_les.create_dataset("u", size_les3d, dtype='float64')
	dset[:,:,:] = np.array(u_in_les)
	
	dset = f_out_les.create_dataset("v", size_les3d, dtype='float64')
	dset[:,:,:] = np.array(v_in_les)
	
	dset = f_out_les.create_dataset("w", size_les3d, dtype='float64')
	dset[:,:,:] = np.array(w_in_les)
	
	f_out_les.close()
	
	f_out_les = h5py.File(fname_out_les, "r")
	eta_out_les = f_out_les["eta"]
	hh_out_les = f_out_les["hh"]
	
	fig2 = plt.figure()
	ax = fig2.add_subplot(3,1,1)
	ax.contourf(eta_out_les)
	ax = fig2.add_subplot(3,1,2)
	ax.contourf(hh_out_les)
	ax = fig2.add_subplot(3,1,3)
	ax.contourf(eta_in_hos)
	
	plt.savefig("projection_fig2.png")
	plt.close()
	
	f_out_les.close()
	

f_in_les.close()
f_in_hos.close()



