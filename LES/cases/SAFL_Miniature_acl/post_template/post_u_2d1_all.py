import sys
import numpy as np
import matplotlib as mpl
#import numpy.fft as npfft
import scipy.fftpack as spfft
from mpi4py import MPI

# A glimpse of this script:
# It is post-process script for numourous snapshots of 2D section of CFD results.
# It has several features:
# 1. use Numpy to store and manipulate the data
# 2. use mpi4py to accelerate computation
# 3. use Matplotlib to plot
# 4. fig1 is timehistory of several locations
# 5. fig2 is mean value profile (1st order moment)
# 6. fig3 is second order moment profile 
# 7. fig4 is quadrant
# 8. fig5 is FFT analysis


# Force matplotlib to not use any Xwindows backend
mpl.use('Agg')
import matplotlib.pyplot as plt

def plotSpectrum(y,Fs):
  """
  Plots a Single-Sided Amplitude Spectrum of y(t)
  """
  n = len(y) # length of the signal
  k = np.arange(n)
  print('n='+str(n)+', Fs='+str(Fs))
  T = float(n)/float(Fs)
  frq = k/T # two sides frequency range
  frq = frq[range(n/2)] # one side frequency range

  #Y = npfft.fft(y)/n # fft computing and normalization
  #Y = Y[range(n/2)]
  
  Y = spfft.fft(y)
  Y = Y[range(n/2)]

  #plot(frq,np.absolute(Y),'r') # plotting the spectrum
  #xlabel('Freq (Hz)')
  #ylabel('|Y(freq)|')
  return Y, frq


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print('Rank assigned as '+str(rank)+'/'+str(size))

isec = 5 
ti_start = 25100 
ti_end = 50000 
ti_interval = 100 
dt = 0.001

ny = 64 
nz = 65

iy0 = 31 

n_ti = (ti_end-ti_start)/ti_interval + 1

temp_mod = np.mod(n_ti, size)
if(rank==0):
  print('n_ti / size = '+str(n_ti)+' / '+str(size)+' = '+str(temp_mod))

if (temp_mod != 0):
  sys.exit('Error: n_ti can not be divided by n_ti')

time = np.zeros(n_ti)
data0 = np.zeros((ny*nz,6))
#data1 = np.zeros((ny*nz,3))
#data2 = np.zeros(nz)

#if rank==0:
#  data_u = np.zeros((n_ti, nz, ny))
#  data_v = np.zeros((n_ti, nz, ny))
#  data_w = np.zeros((n_ti, nz, ny))
#else:
#  data_u = 0
#  data_v = 0
#  data_w = 0

time_local = np.zeros(n_ti/size)
data_u_local = np.zeros((n_ti/size, nz, ny))
data_v_local = np.zeros((n_ti/size, nz, ny))
data_w_local = np.zeros((n_ti/size, nz, ny))

data_u_output = np.zeros((n_ti, 3))  # time, iz
data_u_output_local = np.zeros((n_ti/size, 3))

index = np.zeros((nz, ny), dtype=np.int)

data_velturb_local = np.zeros((n_ti/size, nz, ny, 3))
data_velturb_output_local = np.zeros((n_ti/size, ny, 3))
data_velturb_output = np.zeros((n_ti, ny, 3))
#if rank==0:
#  data_velturb = np.zeros((n_ti, nz, ny, 3))
#data_uuu_mean = np.zeros((nz, 3, 3, 3))

data_uutemp_local = np.zeros((n_ti/size, ny))
data_uusum_local = np.zeros((nz, 3, 3))
data_uusum_global = np.zeros((nz, 3, 3))
  
data_uu_mean = np.zeros((nz, 3, 3))

# data4: mean profile
data_umean = np.zeros(nz)
data_vmean = np.zeros(nz)
data_wmean = np.zeros(nz)

data_usum_global = np.zeros(nz)
data_vsum_global = np.zeros(nz)
data_wsum_global = np.zeros(nz)
  
data_usum_local = np.zeros(nz)
data_vsum_local = np.zeros(nz)
data_wsum_local = np.zeros(nz)

data_umean_2d = np.zeros((nz, ny))
data_vmean_2d = np.zeros((nz, ny))
data_wmean_2d = np.zeros((nz, ny))

data_usum_2d_global = np.zeros((nz, ny))
data_vsum_2d_global = np.zeros((nz, ny))
data_wsum_2d_global = np.zeros((nz, ny))
  
data_usum_2d_local = np.zeros((nz, ny))
data_vsum_2d_local = np.zeros((nz, ny))
data_wsum_2d_local = np.zeros((nz, ny))

print('Memory assigned at cpu '+str(rank)+'/'+str(size))

# data5: timehistory of certain height
iz = 29
iz_array = np.array([28])
#iz_array = np.array([16, 29, 40])
#data5 = np.zeros(n_ti)

for iy in range(ny):
  for i in range(nz):
    index[i, iy] = int(float(iy + i * ny))
  #print(index[i])

#if(rank==0):
  #print(index)

print('------Read files (cpu ' + str(rank) + ')------')

for i in range(n_ti/size):
  ti = ti_start + (n_ti/size*rank + i) * ti_interval
  time_local[i] = ti
  #if ti<=50000:
  #  fname = "../inflow05/output_inlet/inlet_"+"{:010d}".format(ti)+".dat"
  #else:
  #  fname = "../inflow05_2/output_inlet/inlet_"+"{:010d}".format(ti)+".dat"
  #fname = "./POST_U_1D1_"+"{:04d}".format(isec)+"/POST_U_1D1_"+"{:010d}".format(ti)+"_"+"{:04d}".format(isec)+".DAT"
  fname = "POST_U_2D1_"+"{:010d}".format(ti)+"_"+"{:04d}".format(isec)+".DAT"
  data0 = np.genfromtxt(fname, skip_header=2)
  #print(data0)
  if i==0:
    x = data0[0, 0]
    y = data0[index[0,:],1]
    z = data0[index[:,0],2]
  #data1 = data0[:,3:5]
  #data2 = data1[index,0]
  for j in range(nz):
    data_u_local[i, j, :] = data0[index[j,:],3]
    data_v_local[i, j, :] = data0[index[j,:],4]
    data_w_local[i, j, :] = data0[index[j,:],5]

  if(rank==0):
    sys.stdout.write("\rReading files: {:d}/{:d}".format(i+1,n_ti))
    #time.sleep(0.1)
    sys.stdout.flush()

comm.Gather(sendbuf=time_local, recvbuf=time, root=0)
#comm.Gather(sendbuf=data_u_local, recvbuf=data_u, root=0)
#comm.Gather(sendbuf=data_v_local, recvbuf=data_v, root=0)
#comm.Gather(sendbuf=data_w_local, recvbuf=data_w, root=0)


print('\n--------Time averaging at various locations (cpu '+str(rank)+')-------\n')  
for i in range(nz):
  if rank==0:
    sys.stdout.write("\rTime averaging at iz= {:d}/{:d}".format(i+1,nz))
    #time.sleep(0.1)
    sys.stdout.flush()

  data_usum_local[i] = np.sum(data_u_local[:,i,:])
  data_vsum_local[i] = np.sum(data_v_local[:,i,:])
  data_wsum_local[i] = np.sum(data_w_local[:,i,:])
  
  for j in range(ny):
    data_usum_2d_local[i, j] = np.sum(data_u_local[:,i,j])
    data_vsum_2d_local[i, j] = np.sum(data_v_local[:,i,j])
    data_wsum_2d_local[i, j] = np.sum(data_w_local[:,i,j])

comm.Allreduce(sendbuf=data_usum_local, recvbuf=data_usum_global, op=MPI.SUM)
comm.Allreduce(sendbuf=data_vsum_local, recvbuf=data_vsum_global, op=MPI.SUM)
comm.Allreduce(sendbuf=data_wsum_local, recvbuf=data_wsum_global, op=MPI.SUM)

comm.Allreduce(sendbuf=data_usum_2d_local, recvbuf=data_usum_2d_global, op=MPI.SUM)
comm.Allreduce(sendbuf=data_vsum_2d_local, recvbuf=data_vsum_2d_global, op=MPI.SUM)
comm.Allreduce(sendbuf=data_wsum_2d_local, recvbuf=data_wsum_2d_global, op=MPI.SUM)
  
data_umean[:] = data_usum_global[:] / n_ti / ny
data_vmean[:] = data_vsum_global[:] / n_ti / ny
data_wmean[:] = data_wsum_global[:] / n_ti / ny

data_umean_2d[:] = data_usum_2d_global[:] / n_ti 
data_vmean_2d[:] = data_vsum_2d_global[:] / n_ti 
data_wmean_2d[:] = data_wsum_2d_global[:] / n_ti 

#data_wmean[i] = np.average(data_w[:, i, :])

print('\n------Calculate fluctuating components at various locations (cpu '\
  +str(rank)+')-----\n')  

for i in range(nz):
  if rank==0:
    sys.stdout.write("\rCalculate fluctuation at iz= {:d}/{:d}".format(i+1,nz))
    #time.sleep(0.1)
    sys.stdout.flush()
  
  data_velturb_local[:,i,:,0] = data_u_local[:,i,:] - data_umean[i]
  data_velturb_local[:,i,:,1] = data_v_local[:,i,:] - data_vmean[i]
  data_velturb_local[:,i,:,2] = data_w_local[:,i,:] - data_wmean[i]
  for j in range(3):
    for k in range(3):
      data_uutemp_local[:,:] = data_velturb_local[:,i,:,j] * data_velturb_local[:,i,:,k]
      data_uusum_local[i,j,k] = np.sum(data_uutemp_local)

comm.Allreduce(sendbuf=data_uusum_local, recvbuf=data_uusum_global, op=MPI.SUM)
      
data_uu_mean[:,:,:] = data_uusum_global[:,:,:] / n_ti / ny

      #for h in range(3):
      #  temp = data_velturb[:,i,j] * data_velturb[:,i,k] * data_velturb[:,i,h]
      #  data_uuu_mean[i,j,k,h] = np.average(temp)

print('\n')

if(rank==0):
  print('-----Prepare data for plot------')

data_u_output_local = data_u_local[:, iz_array, iy0]
#if(rank==0):
#  print(data_u_output_local)
data_u_output_local=np.ascontiguousarray(data_u_output_local)
#if(rank==0):
#  print(data_u_output_local)
#print(data_u_output_local.flags)
comm.Gather(sendbuf=data_u_output_local, recvbuf=data_u_output, root=0)

data_velturb_output_local = data_velturb_local[:,iz,:,:]
data_velturb_output_local = np.ascontiguousarray(data_velturb_output_local)
comm.Gather(sendbuf=data_velturb_output_local, recvbuf=data_velturb_output, root=0)

if rank == 0 :
  print('----------Plot---------')
## Instantaneous vs z,t
  fig1 = plt.figure()

#ax1_1 = fig1.add_subplot(3,1,1)
#line_u_z1 = ax1_1.plot(time, data_u[:,iz-4])

#ax1_2 = fig1.add_subplot(3,1,2)
#line_u_z2 = ax1_2.plot(time, data_u[:,iz])

#ax1_3 = fig1.add_subplot(3,1,3)
#line_u_z3 = ax1_3.plot(time, data_u[:,iz+5])

  ax1_1 = fig1.add_subplot(1,1,1)
  line_u_z1 = ax1_1.plot(time, data_u_output[:, 0],label='z='+'{:0.3f}'.format(z[iz_array[0]]))
  #line_u_z2 = ax1_1.plot(time, data_u_output[:, 1],label='z='+'{:0.3f}'.format(z[iz_array[1]]))
  #line_u_z3 = ax1_1.plot(time, data_u_output[:, 2],label='z='+'{:0.3f}'.format(z[iz_array[2]]))
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

  um_hub = np.interp(0.125, z, data_umean)
  print("Um_hub="+str(um_hub))

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
  line_u_w = ax4_1.plot(data_velturb_output[:,:,0], data_velturb_output[:,:,2],'b.')
  plt.xlabel("u'")
  plt.ylabel("w'")
  plt.savefig('fig_4_quadrant.png')


## FFT
  fig5 = plt.figure(figsize=(7,7))

# https://jakevdp.github.io/blog/2013/08/28/understanding-the-fft/ 
  #sp = npfft.fft(data_velturb[:,iz,iy0,0])
  #freq = npfft.fftfreq(n_ti) 
  #sp_abs = np.absolute(sp)
  
  Fs = 1000/ti_interval
  sp, freq = plotSpectrum(data_velturb_output[:,iy0,0],Fs)
  
  sp_abs = np.absolute(sp)

  #4:35
    

  ax5_1 = fig5.add_subplot(1,1,1)
  line_ufft_f = ax5_1.plot(freq[1:], sp_abs[1:])
  ax5_1.set_xscale('log')
  ax5_1.set_yscale('log')
  plt.xlabel('freq')
  plt.ylabel('S(f)')
  plt.savefig('fig_5_fft.pdf', format='pdf')
  plt.savefig('fig_5_fft.png')
  np.savez('result_fft',freq, sp_abs)

## Mean vs yz
  fig6 = plt.figure()

  Y, Z = np.meshgrid(y, z, indexing='ij')
  mean_2d = np.transpose(data_umean_2d)
  print(y.shape)
  print(z.shape)
  print(Y.shape)
  print(Z.shape)
  print(data_umean_2d.shape)
  print(mean_2d.shape)
  
  ax6_1 = fig6.add_subplot(1,1,1)
  contour_umean_2d = ax6_1.contourf(Y, Z, mean_2d, 10, cmap=plt.cm.YlOrRd, label='U')
  #ax6_1.legend()
  plt.title('Mean velocity contour')
  plt.xlabel('y')
  plt.ylabel('z')

  #um_hub = np.interp(0.125, z, data_umean)
  #print("Um_hub="+str(um_hub))

  plt.savefig('fig_6_mean_2d.png')

## save mean_2d to tecplot format
  fid = open("umean_2d.dat",'w')
  fid.write('VARIABLES = X, Y, Z, U, V, W\n')
  fid.write('ZONE T="'+str(ti_end * dt)+'" I=1 J='+str(ny)+' K='+str(nz)+' SOLUTIONTIME='\
    +str(ti_end * dt)+'\n')
  for k in range(nz):
    for j in range(ny):
      fid.write(str(x)+' '+str(y[j])+' '+str(z[k])+' '+str(data_umean_2d[k,j])\
        +' '+str(data_vmean_2d[k,j])+' '+str(data_wmean_2d[k,j])+'\n')
  fid.close()

## Save data
  np.savez('result_mean',z,data_umean, data_vmean, data_wmean, data_uu_mean)

  np.savez('result_timehistory', z, data_u_output, iz, iz_array, time, data_velturb_output)
  
  #outfile = 'result_of_post_inlet_5.dat'
  #np.savez(outfile, data_umean, data_vmean, data_wmean,\
  #  data_uu_mean, ti_start, ti_end, ti_interval, ny, nz, iy0, iz, n_ti, time, index)
  
  #np.savez(outfile, data_u, data_v, data_w, data_velturb, data_umean, data_vmean, data_wmean,\
  #  data_uu_mean, ti_start, ti_end, ti_interval, ny, nz, iy0, iz, n_ti, time, index)


#data5 = data_u[:,iz]

#fig1 = plt.figure()
#ax1 = fig1.add_subplot(2,1,1)

#line_zprofile = ax1.plot(data4, z)

#ax2 = fig1.add_subplot(2,1,2)

#line_timehistory = ax2.plot(time, data5)

#plt.savefig('fig_inlet_stat_'+str(ti_start)+'_to_'+str(ti_end)+'.png')




