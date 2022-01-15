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



sx = 7.0  #x 방향 y 방향의 거리.
sy = 7.0
hbar = 500.0
uinfty_o=11.5258407161
uhub_o=9.4546
zhub = 70.0
r = 40.0
D = 2.0*r
PEX = 0.00280499344071

ufric_hi_uw = 0.07495148 # 무슨 파라미터지?

wn_spacing = 2.0*np.pi/(7.0*D*2.0)
wn_swell = 2.0*np.pi/(7.0*D/3.0)

nk = 128 #x방향의 그리드 수 / 3 : 1 side spectrum 의 경우 x 방향 그리드수 의 1/2
# 그리고 spectral method의 dealiasing 2/3 rule 에 의해 나머지 1/3 에 해당하는 높은 주파수 성분은
# 신뢰 할 수 없다. 따라서 전체 x 방향 그리드의 1/3 만 사용.
nz = 129 #z방향의 그리드 수
nk2 = 64 

# we only show part of the spectrum
kindex = np.arange(1, nk2)

casenames = ["les", "hosles", "motion"]
timetags = ["170100to199950", "170100to199950", "170100to199950"]
scaletags = ["original", "scaled"]
#filenames = ["./les/93450to106350/spectrum_mean.dat", "./hosles/147000to166800/spectrum_mean.dat", "./motion/140100to159900/spectrum_mean.dat"]

nvars = [5, 8, 8] #파일의 변수 들의 수. 첫번쨰 케이스는 5개의 변수가 spectrum_mean.dat 파일에 존재
#두번쨰, 세번쨰 케이스는 wave 가 있어, 각각 8개의 변수들이 spectrum_mean.dat 파일에 있다.


###########################################################################
#########다음은 post.f90 에 있는 get_1dspectrum 의 서브 루틴이다 (공부를 위해 여기 일단 추가해본다).#############
sk = 0 # energy sepctral density 즉 energy per one frequency component
    dwk = sqrt((pex*(nx_global/3))**2+(pey*(ny_global/3))**2) / nk # d (wk) 인듯 즉, wave number 의 derivative
   # 즉 pex 는 2pi / LX 이고 여기에 nx_global 을 곱하면 = 2pi / dx 가 된다. 
   # 즉 (pex*(nx_global/3))**2 는 고려 할 수 있는 가장 큰 wave number 이고 이걸 nk로 나누면 dwk 가 됨.
    call fft_for_xy_hos(f,tmp,nx_global,ny_global)
   #fft 를 통해 f (steamwise velocity input) 에 대한 fft 결과 값이 tmp 에 저장.
   #아마 내 생각에는 cross 
    do lp1 = 1, nx_global / 3
       l = lp1 * 2 - 1 # 2를 곱하는 부분은 real part와 imaginary part 이지 않을까 싶음.
       wkx=(l-1)/2*pex # x 방향으로 단위길이당 웨이브가 몇개 있느냐? 숫자가 높아 질수록 파장이 짧아짐. 
   #즉 wkx 는 wave number k in x direction. (0* pex, 1*pex, 2*pex ... nx_global/3 -1 * pex)
   #각 wkx 마다 real part가 있고 imaginary part 가 있을 것임.
       do mp1 = 1, ny_global / 3
          m = mp1 * 2 - 1
          wky=(m-1)/2*pey
          wav=sqrt(wkx**2+wky**2) 
  # 한가지 이해가 안되는 것이. 왜 x y방향을 모두 고려하는건지. z방향에 대해서는 따로 데이터를 저장한다.
  # x, y 방향에 대해서는 단순히 x 방향 단순히 y 방향의 wave number가 아닌, kx ky 의 크기로 
  #k 를 정하고 그 k 가 어디에 해당하는지 다시 아래의 do 문에서 찾은 후에, 에너지를 더 아래에서 더한다.
  # 위의 fft 결과 어떤 값이 생성되는 것인가. tmp 는 kx ky 에 따라 다른 energy spectral density
  # 값을 갖고 있는 것인가 . 예를 들면 temporal fft 에서는 Sxx(f) = |x(f)|^2 으로 나타 낼수 있다.
  # spatial fft 에서는 가령 Sxx(k) = |x(k)|^2 으로 나타나야 할것 같은데 (물론 x(f), x(k) 는 
  # Fourier transform 의 결과 이므로 real part 와 imaginary part 를 함꼐 갖고 있겠지) 
  # 여기서는 spatial fft 결과 tmp는 아마도 Sxy(kx,ky) 의 형태를 지닐 것이다. 그럼 현재 계산 하는 
  # 부분이 cross power spectral density 를 계산 을 하는것인가. 아니면 spatial fft를 통한
  #PSD를 계산 할 때는 원래 이렇게 하는 것인가.
          do n=1,nk
             if(wav >= ((n-0.5)*dwk).and.wav < ((n+0.5)*dwk))then
                exit
             endif
          enddo

          fkn=tmp(l,m)**2+tmp(l,m+1)**2+tmp(l+1,m)**2+tmp(l+1,m+1)**2

          con=4.
          if(l.eq.1.or.m.eq.1) con=2.
          if(l.eq.1.and.m.eq.1) con=1.
          if (n <= nk) sk(n)=sk(n)+con*fkn
       enddo
    enddo

      do n=1,nk
       sk(n)=sk(n)/dwk
    enddo
# 이 일련의 과정을 사실 정확하게 이해 하지 못하겠다. 어떻게 핀은 사용했는지 밑에 좀 보자.
###############################################################################




for iscale in np.arange(2):
	for icase in np.arange(3):
		filename = "./"+casenames[icase]+"/"+timetags[icase]+"/spectrum_mean.dat"
		nvar = nvars[icase]
		data0 = tec.tecplot_reader(filename, [nz, 1, nk, nvar], 2)
#파일네임, 디멘션, 헤더 라인 수. tecplot_reader(filename, [ nz,ny,nx,nvar] ,headerline)
#data 0 를 [nz,nk,nvar]로 reshape. tecplot reader 의 인덱스는 1부터 시작.
# 즉 spectrum_mean.dat 를 nz*1*nk*nvar의 array 로 데이터를 읽어드린다.
#당연히 읽어 드리는 순서는 쓴 순서랑 똑같겠지. 즉. nk 먼저 증가하는 순서대로. 이겠지.
#그래서 읽는 순서를 반대로 적은 거 일수도 있겠네. 즉 nk 가 증가 하는 순으로 먼저 썼기 때문에
#tecplot_reader 에서 읽을 떄는 저 순서대로 읽는 것일 수도 있고. 
#data0 라는 새로운 데이터 파일이 생성된거고. 예를 들어 nz, nk 의 2d 컨투어에 nvar 수의 데이터가 저장 되어있다고 생각하면 된다.

		data0 = data0.reshape([nz, nk, nvar])  #reshape 의 index도 1부터 시작.
		#**이미 data0 에는 spectrum_mean.dat 에 있는 데이터들이 모두 기입 되어 있다.
		#기입 된 방식은 2D 컨투어를 생각하면 편하다.
		##  VARIABLES = X, Z, EK11, EK22, EK33, EKI11, EKI22, EKI33

		z = data0[:,0,1].reshape(nz) 
#저 위에 1 숫자는 data array 를 나타 낸것이므로 저렇게 해야하고 실제 로는 인덱스 0부터 시작하네.
#즉 nk 가 0일 떄 nz 가 0 부터 nz_global -1 까지 그리고 nvar을 z 에 해당하는 1 (두번쨰) 를 
# 읽는다. 읽고 data z 에 저장한다. 그리고 nz 행으로 reshape.
		for i in range(1,nz):
			data0[i, :, 0] = np.arange(nk) * PEX 
#data0 의 0 번쨰 즉 첫번쨰 변수인 x 를 wave number 로 변경해준다. (근데 원래 x 가 wave 번호인데
# 왜 이걸 굳이 다시 변경하지?)
			data0[i, :, 2:] = data0[i, :, 2:]
			if iscale==0:
				data0[i, :, 0] = data0[i, :, 0]
				data0[i, :, 2:] = data0[i, :, 2:] / (ufric_hi_uw**2) 
			elif iscale==1:
				data0[i, :, 0] = data0[i, :, 0] * z[i]
				data0[i, :, 2:] = data0[i, :, 2:] / (ufric_hi_uw**2 * z[i])
##나머지는 tecplot data 파일에서 data extration 을 통해서 그림 그리는게 훨씬 편할듯 하고.
##한가지 물어볼 사항은 scaling 을 어떻게 한건지 반드시 체크.

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
					#x 좌표 가 k[imin1,kindex] kindex 는 1부터 64. y 좌표가 ek. 
					# add some dotted lines to mark the wavenumber of turbine array spacing and swell wavelength and their hamonics.
					ax2.plot([wn_spacing, wn_spacing], [1e-2,1e4], ':')
					# 이 말은 (wn_spacing,1e-1), (wn_spacint,1e4) 이렇게 두 포인트를 :로 이어라.
					ax2.plot([wn_swell, wn_swell], [1e-2,1e4], '--')
					# (wn_sell, 1e-2 ) (wn_swell, 1e4) 이렇게 두 포인트를 -- 로 이어라.
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

