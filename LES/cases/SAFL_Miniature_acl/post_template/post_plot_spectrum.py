import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

from scipy.signal import find_peaks_cwt

import tecplot_io as tec


# for two-column layout, the width is 3.487 inches
fig_width_pt = 300.0  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0/72.27               # Convert pt to inch
golden_mean = (np.sqrt(5)-1.0)/2.0         # Aesthetic ratio
fig_myratio = 1.25
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
          'figure.figsize': fig_size,
          'legend.unicode_minus': False,
          'line.linewidth': 4
}
#          'lines.markerfacecolor': 'none',
#          'scatter.markerfacecolor': 'none'
mpl.rcParams.update(params)


def find_index(_z, _limits):
	_n = len(_z)
	_i_min = 0
	_i_max = _n
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
	

sx = 7.0
sy = 7.0
hbar = 500.0
uinfty_o=11.5258407161
uhub_o=9.4546
zhub = 70.0
r = 40.0
D = 2.0*r
PEX = 0.00280499344071

ufric_hi_uw = 0.07495148

wn_spacing = 2.0*np.pi/(7.0*D*2.0)
wn_swell = 2.0*np.pi/(7.0*D/3.0)

nk = 128
nz = 129
nk2 = 64

# we only show part of the spectrum
kindex = np.arange(1, nk2)

casenames = ["les", "hosles", "motion"]
timetags = ["170100to199950", "170100to199950", "170100to199950"]
scaletags = ["original", "scaled"]
#filenames = ["./les/93450to106350/spectrum_mean.dat", "./hosles/147000to166800/spectrum_mean.dat", "./motion/140100to159900/spectrum_mean.dat"]

nvars = [5, 8, 8]

for iscale in np.arange(2):
	for icase in np.arange(3):
		filename = "./"+casenames[icase]+"/"+timetags[icase]+"/spectrum_mean.dat"
		nvar = nvars[icase]
		data0 = tec.tecplot_reader(filename, [nz, 1, nk, nvar], 2)
		data0 = data0.reshape([nz, nk, nvar])
		##  VARIABLES = X, Z, EK11, EK22, EK33, EKI11, EKI22, EKI33

		z = data0[:,0,1].reshape(nz)

		for i in range(1,nz):
			data0[i, :, 0] = np.arange(nk) * PEX 
			data0[i, :, 2:] = data0[i, :, 2:]
			if iscale==0:
				data0[i, :, 0] = data0[i, :, 0]
				data0[i, :, 2:] = data0[i, :, 2:] / (ufric_hi_uw**2) 
			elif iscale==1:
				data0[i, :, 0] = data0[i, :, 0] * z[i]
				data0[i, :, 2:] = data0[i, :, 2:] / (ufric_hi_uw**2 * z[i])


		k = data0[:,:,0].reshape((nz,nk))
		#print(z)

		zz = np.array((0.03*zhub, zhub - 0.75*D, zhub - 0.5*D, zhub, zhub+0.5*D, 3.0*zhub, 5.0*zhub))
		print("zz = "+str(zz))

		if icase==0:
			eks = np.array((2,3,4))
		else:
			eks = np.array((2,3,4,5,6,7))

		for iek in eks:
			ek = data0[:,:,iek].reshape((nz,nk))
			
			fig = plt.figure()
			ax = fig.add_subplot(1,1,1)

			for izz in range(len(zz)):
				ztemp = zz[izz]
				imin1 = find_index(z, ztemp)
				ztemp = z[imin1]
				#print(data0[imin1, :, [0,2]])
				ax.plot(k[imin1, kindex], ek[imin1, kindex],'-', label=r'$z/H_{hub}='+'{:.2f}'.format(ztemp/zhub)+r'$')
				if iscale==0:
					#indexes = find_peaks_cwt(ek[imin1, kindex], np.arange(1, 10)) 
					#ax.plot(k[imin1, indexes], ek[imin1, indexes], '*')
					if izz==2 and iek == 2:
						for itemp in kindex:
							print(str(itemp)+" "+str(ek[imin1, itemp])+" "+str(k[imin1, itemp])+" "+str(2.0*np.pi/k[imin1, itemp]))
							
					# also we want to plot the spectrum of each elevation individually
					fig2 = plt.figure()
					ax2 = fig2.add_subplot(1,1,1)
					ax2.plot(k[imin1, kindex], ek[imin1, kindex],'-', label=r'$z/H_{hub}='+'{:.2f}'.format(ztemp/zhub)+r'$')
					
					# add some dotted lines to mark the wavenumber of turbine array spacing and swell wavelength and their hamonics.
					ax2.plot([wn_spacing, wn_spacing], [1e-2,1e4], ':')
					ax2.plot([wn_swell, wn_swell], [1e-2,1e4], '--')
					
					plt.xscale('log')
					plt.yscale('log')
					plt.axis('scaled')
					plt.ylim([1e-1,1e3])
					plt.xlabel('$k$')
					plt.ylabel(r'$E_{ij}(z)/(u_*^2)$')
					fig2.tight_layout()
					plt.savefig('spectrum_ek_'+str(iek)+'_'+casenames[icase]+'_'+scaletags[iscale]+'_z_{:03d}'.format(int(np.floor(ztemp/hbar*1000)))+'.png')
					plt.close()
					
					# come back to assemble figure
					plt.figure(fig.number)
					

			if iscale==1:
				xtemp = np.arange(-3.0, 0.0, 0.1)
				xtemp = np.exp(xtemp)
				ytemp = 10**0.8 / xtemp
				ax.plot(xtemp, ytemp, '--')

				xtemp = np.arange(0.1, 4.0, 0.1)
				xtemp = np.exp(xtemp)
				ytemp = 10**0.8 * xtemp**(-5.0/3.0)
				ax.plot(xtemp, ytemp, '--')
				
				ax.plot([1,1],[1e-4,1e2],":")

			#print('k')
			#print(k[imin1, :-1])
			#print('ek11')
			#print(ek11[imin1, :-1])

			plt.xscale('log')
			plt.yscale('log')
			plt.axis('scaled')
			if iscale==0:
				plt.ylim([1e-1,1e3])
				plt.xlabel('$k$')
				plt.ylabel(r'$E_{ij}(z)/(u_*^2)$')
			elif iscale==1 and iek<4:
				plt.ylim([1e-5,1e2])
				#plt.legend()
				plt.xlabel('$kz$')
				plt.ylabel(r'$E_{ij}(z)/(u_*^2 z)$')
			elif iscale==1 and iek>4:
				plt.ylim([1e-4,1e3])
				#plt.legend()
				plt.xlabel('$kz$')
				plt.ylabel(r'$E_{ij}(z)/(u_*^2 z)$')
			

			#if iscale==0 and icase==1 and iek==2:
			#	plt.show()
			
			fig.tight_layout()
			plt.savefig('spectrum_ek_'+str(iek)+'_'+casenames[icase]+'_'+scaletags[iscale]+'.png')
			
			plt.close()

