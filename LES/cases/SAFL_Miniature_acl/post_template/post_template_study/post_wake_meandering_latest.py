import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp

import tecplot_io as tec

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
}
#          'lines.markerfacecolor': 'none',
#          'scatter.markerfacecolor': 'none'
mpl.rcParams.update(params)

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

def plot_bladeline (_ax, _xc, _yc, _r):
	plt.sca(_ax)
	_n = len(_xc)
	for _i in range(_n):
		plt.plot([_xc[_i], _xc[_i]], [_yc[_i]-_r, _yc[_i]+_r], 'k-', linewidth=2)
	return

def get_coord_range (_loc, _D):
	_xc = _loc[0]
	_yc = _loc[1]
	_xran = [_xc+0.2*_D, _xc+6.5*_D]
	_yran = [_yc-2.0*_D, _yc+2.0*_D]
	return _xran, _yran

def get_index_range (_x, _y, _xran, _yran):
	#print "range for x="+str(_x)+" is "+str(_xran)
	#print "range for y="+str(_y)+" is "+str(_yran)
	_ixmin, _ixmax = find_index(_x, _xran)
	_iymin, _iymax = find_index(_y, _yran)
	_ixran = [_ixmin, _ixmax]
	_iyran = [_iymin, _iymax]
	
	return _ixran, _iyran

def collect_timeseries(_ixran, _y_ts, _u, _y2):	
	_iytemp = find_index(_y2, _y_ts)
	_u2 = _u[_ixran[0]:(_ixran[1]+1),_iytemp]
	return _u2	

def identify_wake_core(_ixran, _iyran, _h, _x2, _y2, width=20, _method=2, _hw=6):
	_ix = []
	_iy = []
	_x = []
	_y = []
	#print("ixran="+str(_ixran)+", iyran="+str(_iyran)+", iyran[0]="+str(_iyran[0]))
	
	_tempwidth = (_iyran[1] - _iyran[0])/3
	_tempstart = _iyran[0] + (_iyran[1] - _iyran[0])/2 - _tempwidth/2	
	_tempend = _tempstart + _tempwidth
	
	_utmp = np.zeros(_h.shape[1])
	
	for _i in range(_ixran[0], _ixran[1]+1):
		#_ihmin = np.argmin(_h[_i, _iyran[0]:(_iyran[1]+1)])
		#print("ihmin="+str(_ihmin))
		#_itemp = _iyran[0] + _ihmin
	
		if _method==1:
			# method 1: simply find the global minimum
			_ihmin = np.argmin(_h[_i, _tempstart:_tempend])
		elif _method==2:
			# method 2: firstly get a band average value, then find the location of the minimum
			for _j in range(_iyran[0]+_hw, _iyran[1]-_hw):
				_utmp[_j] = np.average(_h[_i, (_j-_hw):(_j+_hw)])
			_ihmin = np.argmin(_utmp[_tempstart:_tempend])
			
		## I'd like to set the hub as starting point of wake centerline
				
		
		## following code is to improve the continuity of search (only search for the next point in a range of last point.)
		## but it turns out to be easily limited to bottom or top edge of the search frame		
		_itemp = _tempstart + _ihmin
		_tempstart = _itemp - width/2
		
		_tempwidth = width
		
		_tempend = _itemp + width/2
		if _tempstart < (_iyran[0]+_hw):
			_tempstart = _iyran[0] + _hw	
		if _tempend > (_iyran[1]-_hw):
			_tempend = _iyran[1] - _hw
		
		_ix.append(_i)
		_iy.append(_itemp)
		_x.append(_x2[_i])
		_y.append(_y2[_itemp])
		#print(str(_i)+" "+str(_itemp)+" "+str(_x2[_i])+" "+str(_y2[_itemp]))
	return _ix, _iy, _x, _y

def expand_domain(_in, _nepd, mode=0, shift=0):
	_nx = _in.shape[0]
	if len(_in.shape) == 2 :
		_ny = _in.shape[1]
	elif len(_in.shape) > 2 :
		print("Dimension should be 1 or 2.")
		_out = _in
		return _out
		
	_nx2 = _nx + _nepd
	if len(_in.shape) == 1 :
		# 1D case:
		_out = np.empty(_nx2)
		_out[:_nx] = _in[:_nx]
		_out[_nx:_nx2] =  _in[:_nepd]		
	elif len(_in.shape) == 2:
		# 2D case:
		_out = np.empty((_nx2, _ny))
		_out[:_nx, :] = _in[:_nx, :]
		
		if mode == 0:
			_out[_nx:_nx2, :] =  _in[:_nepd, :]	
		elif mode == 1:
			# interpolate/extrapolate
			for _i in range(_ny):				
				_fit = np.polyfit(np.arange(_nx), _in[:_nx, _i] ,1)
				_line = np.poly1d(_fit)
				_out[_nx:_nx2, _i] = _line(np.arange(_nx, _nx2))
				
		if shift == 1 :
			_out2 = _out[_nepd:_nx2, :]
			_out = _out2			
	
	return _out

## http://scipy-cookbook.readthedocs.io/items/SignalSmooth.html
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    #print(len(y))
    #print(len(w))
    #return y
    
    _start = len(w)/2
    return y[(_start):(_start+len(x))]

# These are the "Tableau 20" colors as RGB.  
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120), 
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148), 
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  

# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.  
for i in range(len(tableau20)):  
	r, g, b = tableau20[i]  
	tableau20[i] = (r / 255., g / 255., b / 255.)

sx = 7.0
sy = 7.0
hbar = 500.0
uinfty_o=11.5258407161
uhub_o=9.4546
zhub = 70.0
r = 40.0
D = 2.0*r
PEX = 0.00280499344071
XL = 2.0*np.pi/PEX

L0 = 1.0
U0 = 11.5258407161
T0 = L0/U0
dt = 0.342716489688*T0

ufric = 0.08442747 * U0

nx = 384
ny = 384
nz = 1
nvar = 6

casenames = ["les", "hosles", "motion"]

caseprintnames = ["ONS", "OFX", "OFL"]
#timetags = ["93450to106350", "147000to166800", "170100to174900"]
unametags = ["U", "UI", "UI"]

icase = 0
tis = 103350
tie = 103950
tii = 150
foldername = "./"+casenames[icase]+"/POST_U_2D2_0002/"

y_ts = -2.0/3.0*r

yfftm3 = []
u_ts3 = []

icase = 2

casename_s = 2
casename_e = 2+1
#for icase in range(len(casenames)):
for icase in range(casename_s, casename_e):
	if icase==0:
		tis = 240150
		#tie = 240450+1
		tie = 264900+1
		tii = 150
	elif icase==1:
		tis = 240150
		#tie = 240450+1
		tie = 264900+1
		tii = 150
	elif icase==2:
		tis = 240150
		#tie = 240450+1
		tie = 264900+1
		tii = 150
	else:
		print("icase error.")

	foldername = "./"+casenames[icase]+"/POST_U_2D2_0002/"

	yfftm2 = []
	u_ts2 = []

	turbinefile = "Turbine.inp"
	turbinedata = np.genfromtxt(turbinefile, skip_header=1)
	turbineloc = turbinedata[:,3:5]
	nturbine = len(turbineloc[:,0])

	nepd = 384/8-4 # shift of tecplot grid

	# fig5 is used for plt the overlapped point cloud of all wake velocity minimas.
	#allwakepoints
	fig5, ax5 = plt.subplots()
	xtmp1 = nepd*XL/nx
	xtmp2 = XL + xtmp1
	ax5.plot([xtmp1, xtmp2, xtmp2, xtmp1, xtmp1], [0, 0, XL, XL, 0], 'k-') # draw the domain boundary
	for iturb in range(nturbine):
		xtmp, ytmp = get_coord_range(turbineloc[iturb,:], D)
		ax5.plot([xtmp[0], xtmp[1], xtmp[1], xtmp[0], xtmp[0]], [ytmp[0], ytmp[0], ytmp[1], ytmp[1], ytmp[0]], 'b--') # draw the search frame
	plot_bladeline(ax5, turbineloc[:,0], turbineloc[:,1], r)

	# allmin 
	allmins = []
	for iturb in range(nturbine):
		allmins.append(list())


	for it in range(tis, tie, tii):
		filename = foldername + "POST_"+unametags[icase]+"_2D2_{:010d}".format(it) + "_0002.DAT"
		print(filename)
		data0 = tec.tecplot_reader(filename, [nz, ny, nx, nvar], 2)
		data0 = data0.reshape([ny, nx, nvar])
		##  VARIABLES = X, Y, Z, U, V, W

		X = data0[:,:,0].reshape([ny,nx]).transpose()
		Y = data0[:,:,1].reshape([ny,nx]).transpose()
		U = data0[:,:,3].reshape([ny,nx]).transpose()

		#um = np.average(U)
		#print("um="+str(um))
		um=0.907800731501
		U1 = U / um
		x1 = X[:,0]
		y1 = Y[0,:]
		
		U2 = expand_domain(U1, nepd, shift=1)
		X2 = expand_domain(X, nepd, mode=1, shift=1)
		Y2 = expand_domain(Y, nepd, shift=1)
		x2 = X2[:,0]
		y2 = Y2[0,:]

		
		#print(turbineloc)

		nturbine = len(turbineloc[:,0])
		xran = np.zeros((nturbine,2))
		yran = np.zeros((nturbine,2))
		ixran = np.zeros((nturbine,2), dtype='int')
		iyran = np.zeros((nturbine,2), dtype='int')
		xwake = []
		ywake = []
		ywpo = [] # original location of velocity minima
		ywp = [] # spacial smoothed location of velocity minima
		yfft = []
		yfftm = []
		u_ts = []
		
		for iturb in range(nturbine):
			xran[iturb,:], yran[iturb,:] = get_coord_range(turbineloc[iturb,:], D)
			ixran[iturb,:], iyran[iturb,:] = get_index_range(x2, y2, xran[iturb,:], yran[iturb,:])
			
			y2_ts = (yran[iturb,0] + yran[iturb,1])/2.0 + y_ts
			# u_ts_temp[wake_line_point_index] is velocity at a line along streamline, this line is parallel to rotor centerline with a distance of y_ts
			u_ts_temp = collect_timeseries(ixran[iturb,:], y2_ts, U2, y2)
			# u_ts[turbine_index][wake_line_point_index]
			u_ts.append(u_ts_temp)
			
			ixtemp, iytemp, xtemp, ytemp = identify_wake_core(ixran[iturb,:], iyran[iturb,:], U2, x2, y2)
			xtemp = np.array(xtemp)
			ytemp = np.array(ytemp)
			xwake.append(xtemp)
			ywake.append(ytemp)
			#print("before smooth: len="+str(len(ytemp)))
			#print(ytemp)
			
			ywpo.append(ytemp)
			# allmins[turbine_index][timestep_index][centerlinepoint_index][x,y]
			np_tmp = np.empty((len(ytemp), 2))
			np_tmp[:,0] = xtemp
			np_tmp[:,1] = ytemp
			allmins[iturb].append(np_tmp)
			
			ytemp = smooth(ytemp)
			#print("after smooth: len="+str(len(ytemp)))
			#print(ytemp)
			ywp.append(ytemp)
			#print("iturb="+str(iturb))
			#print(np.array(ixtemp))
			#print(np.array(iytemp))
			#print(np.vstack((np.array(ixtemp), np.array(iytemp))))
			
			p_temp = np.abs(np.fft.fft(ytemp-np.average(ytemp)))
			yfft.append(p_temp)	
			if iturb==0:
				yfftm = p_temp
			else:			
				yfftm = yfftm + p_temp
			
		# u_ts2[timestep_index][turbine_index][centerlinepoint_index]
		u_ts2.append(u_ts)
		
		xfft = 2.0*np.pi/(x2[1]-x2[0]) * np.arange(len(yfftm))
		yfftm = yfftm / nturbine
		if it==tis:
			yfftm2 = yfftm
		else:
			yfftm2 = yfftm2 + yfftm

		fig, ax = plt.subplots()
		#levels = np.arange(0.5, 1.5, 0.01) 
		#ax1 = plt.contourf(X, Y, U2, level=levels, cmap=plt.get_cmap(name ="Oranges_r"), label='U')

		normalize = mpl.colors.Normalize(vmin=U2.min(), vmax=U2.max())
		#cmap=plt.get_cmap("Oranges_r")
		cmap=plt.get_cmap("coolwarm_r")

		color = U2
		#color=cmap(normalize(U2[i]))

		## 1D color input
		color = U2
		#color = cmap((color-color.min())/(color.max()-color.min()))
		#color = np.array(color)

		## 2D color input
		#color = U2
		#cnp = np.array(color)
		#color = np.empty([color.shape[0],color.shape[1],4])
		#for i in range(color.shape[0]):
		#	for j in range(color.shape[1]):
		#		color[i,j,:] = cmap((cnp[i,j]-cnp[:,:].min())/(cnp[:,:].max()-cnp[:,:].min()))

		#ax1 = plt.contourf(X, Y, U2, cmap=color, label='U')

		#ax1 = plt.contourf(X2, Y2, U2, cmap=cmap, normal=normalize, label='U')
		#cbar = fig.colorbar(ax1)
		
		uc = (0.25*(U2[:-1,:-1] + U2[:-1,1:] + U2[1:,:-1] + U2[1:,1:]))
		ax1 = plt.pcolor(X2,Y2,uc, cmap=cmap, vmin=0.3, vmax=1.6)
		cbar = plt.colorbar(ax1)
		plt.text(1200, -150, "$tu_*/H="+"{:.3f}".format(it * dt * ufric / hbar)+"$")

		plot_bladeline(ax, turbineloc[:,0], turbineloc[:,1], r)

		for iturb in range(nturbine):
			xran2 = x2[ixran[iturb,:]]
			yran2 = y2[iyran[iturb,:]]
			plt.plot([xran2[0],xran2[1],xran2[1],xran2[0],xran2[0]],[yran2[0],yran2[0],yran2[1],yran2[1],yran2[0]],'b--')
			plt.plot(xwake[iturb], ywake[iturb], 'o', color=tableau20[2])
			plt.plot(xwake[iturb], ywp[iturb], 'k-', linewidth=3, alpha=0.5)
			
			plt.figure(fig5.number)
			ax5.plot(xwake[iturb], ywake[iturb], '.', color=tableau20[1], alpha=0.5)
			
			plt.figure(fig.number)

		#plt.show()
		plt.xlim([x2.min(), x2.max()])
		plt.ylim([y2.min(), y2.max()])
		#plt.axis('equal')
		plt.savefig(foldername+"wake_{:010d}".format(it)+".png")
		plt.close()

		#fig2 = plt.figure()
		#for iturb in range(nturbine):
		#	ax2 = fig2.add_subplot(4,4,iturb)
		#	plt.plot(xfft, yfft[iturb])
		#	plt.xscale('log')
		#	plt.yscale('log')
		#	#plt.axis('scaled')
		#plt.show()
		#plt.close()

	plt.figure(fig5.number)
	plt.savefig(foldername+"allwake.png")
	plt.close()

	yfftm2 = yfftm2 / (len(range(tis, tie, tii)))
	yfftm3.append(yfftm2)
	fig3 = plt.figure()
	ax3 = fig3.add_subplot(1,1,1)
	plt.plot(xfft, yfftm2)
	plt.xscale('log')
	plt.yscale('log')
	plt.ylim([1,1e4])
	plt.xlabel('$k$')
	plt.ylabel(r'$S_y$')
	plt.savefig(foldername+"wake_spectrum.png")
	plt.close()

	np.savetxt(foldername+"wave_spectrum.dat", np.vstack((xfft, yfftm2)))

	#u_ts3.append(u_ts2)
	lentemp = np.empty(nturbine, dtype='int')
	for i in range(nturbine):
		temp = u_ts2[0][i]
		lentemp[i] = len(temp)
	xlen = np.min(lentemp)
	tlen = len(range(tis,tie,tii))

	u_ts_o = np.empty((nturbine, xlen, tlen))
	for i in range(nturbine):
		for j in range(tlen):
			u_ts_o[i, 0:xlen, j] = u_ts2[j][i][0:xlen]
			
	np.savez(foldername+"u_ts.npz", u_ts_o)
	n_ts_o = []

	## save the wake_velocity_minima locations:
	wkp_all = np.empty((nturbine, tlen, xlen, 2))
	for i in range(nturbine):
		for j in range(tlen):
			wkp_all[i, j, :, :] = allmins[i][j][0:xlen,:]
	np.savez(foldername+"wake_minima_points.npz", wkp_all=wkp_all)
	wkp_all = []
		
	
			
	
#fig4 = plt.figure()
#ax4 = fig4.add_subplot(1,1,1)
#for i in range(len(casenames)):
#	plt.plot(xfft, yfftm3[i], '-', label=caseprintnames[i])
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim([1,1e4])
#plt.xlabel('$k$')
#plt.ylabel(r'$S_y$')
#plt.legend()
#plt.savefig("wake_spectrum_compare.png")
#plt.close()
