import sys
import numpy as np
import matplotlib as mpl
import numpy.fft as npfft
# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

ti_start = 35000 
ti_end = 74000 
ti_interval = 1 

ny = 64 
nz = 65

iy0 = 40 

n_ti = (ti_end-ti_start)/ti_interval + 1

time = np.zeros(n_ti)
data0 = np.zeros((ny*nz,6))
data1 = np.zeros((ny*nz,3))
data2 = np.zeros(nz)
data_u = np.zeros((n_ti, nz, ny))
data_v = np.zeros((n_ti, nz, ny))
data_w = np.zeros((n_ti, nz, ny))
data_velturb = np.zeros((n_ti, nz, ny, 3))
#data_uuu_mean = np.zeros((nz, 3, 3, 3))
data_uu_mean = np.zeros((nz, 3, 3))

index = np.zeros((nz, ny), dtype=np.int)

# data4: mean profile
data_umean = np.zeros(nz)
data_vmean = np.zeros(nz)
data_wmean = np.zeros(nz)

# data5: timehistory of certain height
iz = 15 
#data5 = np.zeros(n_ti)

for iy in range(ny):
  for i in range(nz):
    index[i, iy] = int(float(iy + i * ny))
  #print(index[i])

print(index)

print('------Read files------')

for i in range(n_ti):
  ti = ti_start + i * ti_interval
  time[i] = ti
  if ti<=50000:
    fname = "../inflow05/output_inlet/inlet_"+"{:010d}".format(ti)+".dat"
  else:
    fname = "output_inlet/inlet_"+"{:010d}".format(ti)+".dat"
  data0 = np.genfromtxt(fname, skip_header=2)
  #print(data0)
  if i==0:
    z = data0[index[:,0],2]
  data1 = data0[:,3:5]
  #data2 = data1[index,0]
  for j in range(nz):
    data_u[i, j, :] = data0[index[j,:],3]
    data_v[i, j, :] = data0[index[j,:],4]
    data_w[i, j, :] = data0[index[j,:],5]

  sys.stdout.write("\rReading files: {:d}/{:d}".format(i,n_ti))
  #time.sleep(0.1)
  sys.stdout.flush()


  #print(data3)

#np.savetxt('inlet_data_u_'+str(ti_start)+'_to_'+str(ti_end)+'.dat',data_u)
#np.savetxt('inlet_data_v_'+str(ti_start)+'_to_'+str(ti_end)+'.dat',data_v)
#np.savetxt('inlet_data_w_'+str(ti_start)+'_to_'+str(ti_end)+'.dat',data_w)

print('--------Process data at various locations-------')  
for i in range(nz):
  sys.stdout.write("\rProcessing data: {:d}/{:d}".format(i,nz))
  #time.sleep(0.1)
  sys.stdout.flush()

  data_umean[i] = np.average(data_u[:, i, :])
  data_vmean[i] = np.average(data_v[:, i, :])
  data_wmean[i] = np.average(data_w[:, i, :])
  data_velturb[:,i,:,0] = data_u[:,i,:] - data_umean[i]
  data_velturb[:,i,:,1] = data_v[:,i,:] - data_vmean[i]
  data_velturb[:,i,:,2] = data_w[:,i,:] - data_wmean[i]
  for j in range(3):
    for k in range(3):
      temp = data_velturb[:,i,:,j] * data_velturb[:,i,:,k]
      data_uu_mean[i,j,k] = np.average(temp)

      #for h in range(3):
      #  temp = data_velturb[:,i,j] * data_velturb[:,i,k] * data_velturb[:,i,h]
      #  data_uuu_mean[i,j,k,h] = np.average(temp)



## Instantaneous vs z,t
fig1 = plt.figure()

#ax1_1 = fig1.add_subplot(3,1,1)
#line_u_z1 = ax1_1.plot(time, data_u[:,iz-4])

#ax1_2 = fig1.add_subplot(3,1,2)
#line_u_z2 = ax1_2.plot(time, data_u[:,iz])

#ax1_3 = fig1.add_subplot(3,1,3)
#line_u_z3 = ax1_3.plot(time, data_u[:,iz+5])

ax1_1 = fig1.add_subplot(1,1,1)
line_u_z1 = ax1_1.plot(time, data_u[:, iz-4, iy0],label='z='+'{:0.3f}'.format(z[iz-4]))
line_u_z2 = ax1_1.plot(time, data_u[:,   iz, iy0],label='z='+'{:0.3f}'.format(z[iz]))
line_u_z3 = ax1_1.plot(time, data_u[:, iz+5, iy0],label='z='+'{:0.3f}'.format(z[iz+5]))
ax1_1.legend()
plt.title('Streamwise velocity at various heights')
plt.xlabel('timestep')
plt.ylabel('u')

plt.savefig('fig_1_instantaneous_various_z.png')


## Mean vs z
fig2 = plt.figure()

ax2_1 = fig2.add_subplot(1,1,1)
line_umean_z = ax2_1.plot(data_umean, z, label='U')
line_vmean_z = ax2_1.plot(data_vmean, z, label='V')
line_wmean_z = ax2_1.plot(data_wmean, z, label='W')
ax2_1.legend()
plt.title('Mean velocity vs. height')
plt.xlabel('Mean Velocity')
plt.ylabel('z')

plt.savefig('fig_2_mean_z.png')


## u'_i * u'_j vs. z
fig3 = plt.figure()

ax3_1 = fig3.add_subplot(1,1,1)
#plt.rc('text', usetex=True)
#plt.rc('font', family='Arial')

line_uu_z = ax3_1.plot(data_uu_mean[:,0,0], z, label=r"$\overline{u^\prime u^\prime}$")
line_vv_z = ax3_1.plot(data_uu_mean[:,1,1], z, label=r"$\overline{v^\prime v^\prime}$")
line_ww_z = ax3_1.plot(data_uu_mean[:,2,2], z, label=r"$\overline{w^\prime w^\prime}$")
line_uw_z = ax3_1.plot(-data_uu_mean[:,0,2], z, label=r"$\overline{-u^\prime w^\prime}$")

#line_uu_z = ax3_1.plot(data_uu_mean[:,0,0], z, label=r"$u'u'$")
#line_vv_z = ax3_1.plot(data_uu_mean[:,1,1], z, label=r"$v'v'$")
#line_ww_z = ax3_1.plot(data_uu_mean[:,2,2], z, label=r"$w'w'$")
#line_uw_z = ax3_1.plot(data_uu_mean[:,0,2], z, label=r"$u'w'$")

ax3_1.legend()
plt.title('Second order velocity moments vs. height')
plt.xlabel('2nd order moments')
plt.ylabel('z')

plt.savefig('fig_3_2ndmoment_z.png')


## quadrant
fig4 = plt.figure()
ax4_1 = fig4.add_subplot(1,1,1)
line_u_w = ax4_1.plot(data_velturb[:,iz,:,0], data_velturb[:,iz,:,2],'b.')
plt.xlabel("u'")
plt.ylabel("w'")
plt.savefig('fig_4_quadrant.png')


## FFT
fig5 = plt.figure()

sp = npfft.fft(data_velturb[:,iz,iy0,0])
freq = npfft.fftfreq(n_ti)

sp_abs = np.absolute(sp)

ax5_1 = fig5.add_subplot(1,1,1)
line_ufft_f = ax5_1.plot(freq, sp_abs)
plt.xlabel('freq')
plt.ylabel('FFT')
plt.savefig('fig_5_fft.png')


#data5 = data_u[:,iz]

#fig1 = plt.figure()
#ax1 = fig1.add_subplot(2,1,1)

#line_zprofile = ax1.plot(data4, z)

#ax2 = fig1.add_subplot(2,1,2)

#line_timehistory = ax2.plot(time, data5)

#plt.savefig('fig_inlet_stat_'+str(ti_start)+'_to_'+str(ti_end)+'.png')




