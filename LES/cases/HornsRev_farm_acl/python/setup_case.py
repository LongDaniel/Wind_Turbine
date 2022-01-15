import numpy as np

# This script aims to automatically generate input files for LES_HOS_TURBINE simulation

## basic parameters, some values might be overwriten in following layout setups
rotor_model = 3 
nacelle_model = 1

NP1 = 16 
NP2 = 16
NX = 384 
NY = 384
NZ = 129

NSEG = 3
# 12 82 _ 0.008 _ _ 6.0 5.0 2.0
NZ1 = 6 
NZ2 = 30 
NZ3 = NZ - NZ1 - NZ2 
L1 = 0.05
RAT1 = 1.1 
RAT2 = 1.4 
RAT3 = 1.522 

nstep_oneperiod = 125.0
nperiod = 40.0
DT0 = 0.001


NOUTD = 150 
NOUTC = 150

CLBETA = 1.0
CLGAMA = 1.0

ARM = 0.08

Shen1_AL = 1 
Shen1_AL_tipcorrection = 1.0
Shen1_AL_tipcorrectionratio_Fa = 1.0

InletRelaxation = 0 
IR_start = 4.15
IR_end = 4.2
IR_cx = 0.1
IR_cy = 0.1
IR_cz = 0.1
IR_NSmooth=0
IR_dataindex = 1
# IR_scale = 0.183730219

c_bforce = 0.75

adjust_CT_ACL = 0 
adjust_CT_nacelle = 2 
adjust_CP_ACL = 0
CT0_ACL = 0.75
CP0_ACL = 0.20
CT0_nacelle = 0.5
ratio_CT0 = np.ones(3)
WE = 5.0e-4

mode_inflow = 0 

# default value:
iwavy = 1 
AKA = 0.0
NSWAVEX = 4
USTAR = 0.029604
U10 = 8.0
USS = 0.44
GAMMA = 3.3
F = 15000.0
G = 9.81
PHI = 90
IST = 0
periodp_scl = 10.0
NXMOD = NX
NYMOD = NY

layout = 7 
# some preset layouts:
# 0, 
# 1, single Vestas V80 
# 2, single Clipper
# 3, single GWS/EP-6030x3
# 5, Horns Rev 1 wind farm, 10x7 Vestas V80 array
# 6, Horns Rev 1 wind farm, 3x3 array, real wave
# 7, Horns Rev 1 Wind Farm, 3x3 array, ideal swell

if layout==1 :
  d_o = 80.0 # rotor diameter
  rmin_o = 0.07507508*40.0 # blade root distance to axis
  hhub_o = 70.0 # turbine hub height
  uhub_o = 8.0 # wind velocity at hub height (undistrubed or at rotor plane?)
  omega_o = 1.69047 # turbine angular velocity in rad/s
  tsr_o = 8.45235 # Tip speed ratio
  pitch_o = -1.02 # pitch in degree
  nfoiltype = 1
  
  spacing_x_o = 7.0 * d_o
  spacing_incl = 7.0 / 180.0 * np.pi
  spacing_y_o = spacing_x_o * np.cos(spacing_incl)

  nturbine_x = 1
  nturbine_y = 1
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 500.0 # boundary layer height of atmosphere
  utao_air_o = 0.442 # boundary layer friction velocity
  z0_air_o = 0.05 # boundary layer friction height
  kappa = 0.4 # von Karmon number
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = spacing_x_o
  domain_y_o = spacing_y_o
  domain_z_o = delta_air_o
if layout==5 :
  d_o = 80.0 # rotor diameter
  rmin_o = 0.033694/0.448799*40.0 # blade root distance to axis
  hhub_o = 70.0 # turbine hub height
  uhub_o = 8.0 # wind velocity at hub height (undistrubed or at rotor plane?)
  omega_o = 1.69047 # turbine angular velocity in rad/s
  tsr_o = 8.45235 # Tip speed ratio
  pitch_o = -1.02 # pitch in degree
  nfoiltype = 1

  inflowdir = 0.0 / 180.0 * np.pi
  
  spacing_x_o = 7.0 * d_o * np.cos(inflowdir)
  spacing_incl = 7.0 / 180.0 * np.pi
  spacing_y_o = 7.0 * d_o * np.cos(spacing_incl-inflowdir)

  nturbine_x = 3 
  nturbine_y = 3
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 500.0 # boundary layer height of atmosphere
  utao_air_o = 0.442 # boundary layer friction velocity
  z0_air_o = 0.05 # boundary layer friction height
  kappa = 0.4 # von Karmon number
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = 2.0 * nturbine_x * spacing_x_o
  domain_y_o = 1.5 * nturbine_y * spacing_y_o
  domain_z_o = delta_air_o
# wind farm
if layout==6 :
  ihos = 1
  iwavy = 3
  
  d_o = 80.0 # rotor diameter
  rmin_o = 0.033694/0.448799*40.0 # blade root distance to axis
  hhub_o = 70.0 # turbine hub height
  uhub_o = 9.4546 # wind velocity at hub height (undistrubed or at rotor plane?)
  omega_o = 1.69047 # turbine angular velocity in rad/s
  tsr_o = 8.45235 # Tip speed ratio
  pitch_o = -1.02 # pitch in degree
  nfoiltype = 1

  inflowdir = 0.0 / 180.0 * np.pi
  
  spacing_x_o = 7.0 * d_o * np.cos(inflowdir)
  spacing_incl = 7.0 / 180.0 * np.pi
  spacing_y_o = 7.0 * d_o * np.cos(spacing_incl-inflowdir)

  nturbine_x = 3 
  nturbine_y = 3
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 500.0 # boundary layer height of atmosphere
  utao_air_o = 0.43190 # boundary layer friction velocity
  z0_air_o = 0.00885475 # boundary layer friction height
  kappa = 0.41 # von Karmon number
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = 1.0 * nturbine_x * spacing_x_o
  domain_y_o = 1.0 * nturbine_y * spacing_y_o
  domain_z_o = delta_air_o

  CT0_ACL = 0.60

  NXMOD = 256 
  NYMOD = 256
  NPW = 3 
  NCPU = NP1 * NP2
  NSWAVEX = 24

  #Jonswap
  OMEGA_P = 0.9322
  GAMMA = 3.3
  F = 170818.0
  G = 9.81
  PHI = 90

  # surface tension
  IST = 0
  ICAP = 0
  SIGMA = 0.0728

  ISWELL = 0
  AKA = 0.645 * 2.0 * np.pi / (domain_x_o/NSWAVEX)
  NSWELLX = 8
  NSWELLY = 8

  NTP = 128 
  DT_HOS0 = -0.2
  U0_HOS = uinfty_o
  L0_HOS = 1.0
  G0_HOS = 9.81
  FR2 = U0_HOS * U0_HOS / G0_HOS / L0_HOS
  USS = 0.45
  U10 = 4.17
  USTAR = USS / U0_HOS 

  NTIME_HOS = 64000

  NFHOS = 128
  NTHOS = 128
  
  omegap_o = np.sqrt(NSWAVEX*2.0*np.pi/(domain_x_o)*G0_HOS) 
  periodp_o = np.pi * 2.0 / omegap_o
  omegap_scl = np.sqrt(NSWAVEX*2.0*np.pi/(domain_x_o/L0_HOS)/FR2) 
  periodp_scl = np.pi * 2.0 / omegap_scl
  print("Original: Omega_p = "+str(omegap_o)+", Period_p = "+str(periodp_o))
  print("Scaled: Omega_p = "+str(omegap_scl)+", Period_p = "+str(periodp_scl))

if layout==7 :
  ihos = 1
  iwavy = 7 
  
  d_o = 80.0 # rotor diameter
  rmin_o = 0.033694/0.448799*40.0 # blade root distance to axis
  hhub_o = 70.0 # turbine hub height
  uhub_o = 9.4546 # wind velocity at hub height (undistrubed or at rotor plane?)
  omega_o = 1.69047 # turbine angular velocity in rad/s
  tsr_o = 8.45235 # Tip speed ratio
  pitch_o = -1.02 # pitch in degree
  nfoiltype = 1

  inflowdir = 0.0 / 180.0 * np.pi
  
  spacing_x_o = 7.0 * d_o * np.cos(inflowdir)
  spacing_incl = 0.0 / 180.0 * np.pi
  spacing_y_o = 7.0 * d_o * np.cos(spacing_incl-inflowdir)

  nturbine_x = 4 
  nturbine_y = 4
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 500.0 # boundary layer height of atmosphere
  utao_air_o = 0.43190 # boundary layer friction velocity
  z0_air_o = 0.00885475 # boundary layer friction height
  kappa = 0.41 # von Karmon number
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = 1.0 * nturbine_x * spacing_x_o
  domain_y_o = 1.0 * nturbine_y * spacing_y_o
  domain_z_o = delta_air_o

  CT0_ACL = 0.60

  NXMOD = 768 
  NYMOD = 768
  NPW = 3 
  NCPU = NP1 * NP2
  NSWAVEX = 12 

  #Jonswap
  OMEGA_P = 0.9322
  GAMMA = 3.3
  F = 170818.0
  G = 9.81
  PHI = 90

  # surface tension
  IST = 0
  ICAP = 0
  SIGMA = 0.0728

  ISWELL =1 
  #AKA = 0.645 * 2.0 * np.pi / (domain_x_o/NSWAVEX)
  AKA = 0.1
  NSWELLX = 12
  NSWELLY = 1

  NTP = 128 
  DT_HOS0 = -0.2
  U0_HOS = uinfty_o
  L0_HOS = 1.0
  G0_HOS = 9.81
  FR2 = U0_HOS * U0_HOS / G0_HOS / L0_HOS
  USS = 0.45
  U10 = 4.17
  USTAR = USS / U0_HOS 

  NTIME_HOS = 64000

  NFHOS = 128
  NTHOS = 128
  
  omegap_o = np.sqrt(NSWAVEX*2.0*np.pi/(domain_x_o)*G0_HOS) 
  periodp_o = np.pi * 2.0 / omegap_o
  omegap_scl = np.sqrt(NSWAVEX*2.0*np.pi/(domain_x_o/L0_HOS)/FR2) 
  periodp_scl = np.pi * 2.0 / omegap_scl
  print("Original: Omega_p = "+str(omegap_o)+", Period_p = "+str(periodp_o))
  print("Scaled: Omega_p = "+str(omegap_scl)+", Period_p = "+str(periodp_scl))

  fsitype = 1 
  fsidof = np.array((1,0,1,0,1,0))
  fsirefpos_o = np.zeros((nturbine_x, nturbine_y, 6)) # it should be updated later
  fsiamp_o = np.array((0.74272, 0.0, 0.14854, 0.0, 4.0*np.pi/180.0, 0.0))
  fsiomega_o = omegap_o * np.ones(6)
  fsiphase0 = np.pi * np.array((0.0, 0.0, 1.0, 0.0, 1.0, 0.0)) 

  
   
elif layout==2 :
  d_o = 1.0 # rotor diameter
  rmin_o = 2.8/96.0 # blade root distance to axis
  hhub_o = 2.0 # turbine hub height
  uhub_o = 1.0 # wind velocity at hub height (undistrubed or at rotor plane?)
  omega_o = 10.0 # turbine angular velocity in rad/s
  tsr_o = 5.0 # Tip speed ratio
  pitch_o = 0.0 # pitch in degree
  nfoiltype = 5
  
  spacing_x_o = 7.0 * d_o
  spacing_incl = 0.0 / 180.0 * np.pi
  spacing_y_o = spacing_x_o * np.cos(spacing_incl)

  nturbine_x = 1
  nturbine_y = 1
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 4.0 # boundary layer height of atmosphere
  z0_air_o = 0.0001 * delta_air_o # boundary layer friction height
  kappa = 0.4 # von Karmon number
  utao_air_o1 = uhub_o*kappa/np.log(hhub_o/z0_air_o) # boundary layer friction velocity
  #utao_air_o2 = 0.102
  utao_air_o = utao_air_o1
  #print('utao_air_o: '+str(utao_air_o1)+' '+str(utao_air_o2))
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = spacing_x_o
  domain_y_o = delta_air_o 
  domain_z_o = delta_air_o
elif layout==3 :
  d_o = 0.15 # rotor diameter
  rmin_o = 0.01 # blade root distance to axis
  hhub_o = 0.125 # turbine hub height
  uhub_o = 2.2 # wind velocity at hub height (undistrubed or at rotor plane?)
  tsr_o = 4.4 # Tip speed ratio
  omega_o = tsr_o * uhub_o / (d_o / 2.0)  # turbine angular velocity in rad/s
  pitch_o = 0.0 # pitch in degree
  nfoiltype = 1
  
  spacing_x_o = 4.32
  spacing_incl = 0.0 / 180.0 * np.pi
  spacing_y_o = 0.72

  nturbine_x = 1
  nturbine_y = 1
  nturbine = nturbine_x * nturbine_y
  
  rho_air = 1.205
  nu_air = 1.511e-5
  
  delta_air_o = 0.46 # boundary layer height of atmosphere
  z0_air_o = 0.000018 # boundary layer friction height
  kappa = 0.41 # von Karmon number
  utao_air_o1 = uhub_o*kappa/np.log(hhub_o/z0_air_o) # boundary layer friction velocity
  utao_air_o2 = 0.102
  utao_air_o = utao_air_o1
  print('utao_air_o: '+str(utao_air_o1)+' '+str(utao_air_o2))
  
  uinfty_o = utao_air_o / kappa *np.log(delta_air_o / z0_air_o) # U_infty, guessed through log law
  re_tao_o = utao_air_o * delta_air_o /nu_air # friction Reynolds number
  re_o = uinfty_o * delta_air_o / nu_air # Reynolds number

  domain_x_o = spacing_x_o
  domain_y_o = spacing_y_o 
  domain_z_o = delta_air_o
  
if mode_inflow == 1:
  #domain_x_o = 2.0*domain_z_o
  rotor_model =  -100
  NOUTC = 1
  NOUTD = int(nstep_oneperiod * 1)
  NX = int(2.0*(NZ-1))
  nperiod = 200


    

u0 = uinfty_o
#l0 = spacing_x_o / (2.0 * np.pi)
#u0 = 1.0
l0 = 1.0
g0 = 9.81
RE = u0*l0/nu_air
domain_x_scl = domain_x_o / l0
domain_y_scl = domain_y_o / l0
domain_z_scl = domain_z_o / l0

print("Uinfty_o="+str(uinfty_o)+', uhub_o='+str(uhub_o))

NTIME = int(nperiod * nstep_oneperiod)
hhub_scl = hhub_o / l0
d_scl = d_o / l0
radius_scl = d_o / 2.0 / l0
print("Diameter scale: from "+str(d_o)+" to "+str(2.0*radius_scl)+" (L0="+str(l0)+")")

dx_scl = domain_x_scl / NX
dy_scl = domain_y_scl / NY
dz_scl = domain_z_scl / NZ
print("Choose of nx, dx:")
tmp = "Lx="+str(domain_x_scl)+", Nx="+str(NX)+", dx/D="+str(dx_scl / d_scl)
print(tmp)
tmp = "Ly="+str(domain_y_scl)+", Ny="+str(NY)+", dy/D="+str(dy_scl / d_scl)
print(tmp)
tmp = "Lz="+str(domain_z_scl)+", Nz="+str(NZ)+", dz/D="+str(dz_scl / d_scl)
print(tmp)
tmp = "LZ1="+str(hhub_scl+radius_scl)+", Nz="+str(NZ1)+", dz/D="+str((hhub_scl+radius_scl)/NZ1/ d_scl)
print(tmp)
print("NZ1 is good around " + str((hhub_scl+radius_scl)/dy_scl))

PEX = 2.0 * np.pi / domain_x_scl
PEY = 2.0 * np.pi / domain_y_scl
ZL = 1.0
HBAR = domain_z_scl
Z0 = z0_air_o / l0
RESTOP = 0.0
RESBOT = RE / (1.0 / kappa * np.log(HBAR / Z0))
FR2 = u0*u0/g0/l0
omega_scl = omega_o / (u0 / l0)
period_scl = 2.0 * np.pi / omega_scl

# scaling for HOS parameters
if ihos==1:
  BOND = 1000.0/(l0*l0)/SIGMA


# Velocity profile: u(z)=us/kappa*log(z/z0), where us=1/(2.5*log(hbar/z0))
ztemp = np.zeros(20)
utemp = np.zeros(20)
for i in range(20):
  ztemp[i]=(i+1.0)/20.0*HBAR
  utemp[i]=utao_air_o*2.5*np.log(ztemp[i]/z0_air_o)

z_and_u = np.column_stack((ztemp,utemp))
print('z vs. u')
print(z_and_u)

L12 = (hhub_scl + radius_scl)/HBAR

grid_res = np.sqrt(dy_scl*dz_scl)
DT = period_scl / nstep_oneperiod
if layout==3:
  DT=DT0

dt_x = dx_scl / 1.0
dt_yz = grid_res / (omega_scl * radius_scl)
print("\nchoose of dt: preset dt="+ str(DT))
print("dx="+str(dx_scl)+", ux=1.0, dt_x = "+str(dt_x))
print("dy="+str(grid_res)+", urot="+str(omega_scl*radius_scl)+", dt_yz="+str(dt_yz))

print("Runtime="+str(DT)+"x"+str(NTIME)+"="+str(DT*NTIME)+"="+str(DT*NTIME/period_scl)+"T"\
  +" (T="+str(period_scl)+")")

## fort.11
fid = open("fort.11.template",'w')
fid.write("0 # ISTART\n")
fid.write("0 # ISCALAR\n")
fid.write(str(rotor_model)+" # ITURBINE\n")
fid.write("0 # IWCONTROL\n")
fid.write(str(NP1)+" "+str(NP2)+" # NP1, NP2\n")
fid.write(str(PEX)+" "+str(PEY)+" # PEX, PEY\n")
fid.write(str(ZL)+" "+str(HBAR)+" # ZL, HBAR\n")
fid.write(str(Z0)+" # Z0\n")
fid.write(str(NX)+" "+str(NY)+" "+str(NZ)+" # NX, NY, NZ\n")
fid.write(str(NXMOD)+" "+str(NYMOD)+" # NXS, NYS\n")
fid.write(str(DT)+" "+str(NTIME)+" # DT, NTIME\n")
fid.write("20 1.E-8 # ITMAX, ERLIM\n")
fid.write(str(RESBOT)+" "+str(RESTOP)+" # RESBOT, RESTOP\n")
fid.write(str(FR2)+" "+str(WE)+" # FR2, WE\n")
fid.write(str(NOUTD)+" "+str(NOUTC)+" # NOUTD, NOUTC\n")
fid.write(str(CLBETA)+" "+str(CLGAMA)+" # CLBETA, CLGAMA\n")
fid.write(str(ARM)+" # ARM\n")
fid.write("0.0 3 2 # HKA, NWAVE, CRAT\n")
fid.write("1.0 0.1 # TIMEWAVY, TCOEF\n")
fid.write(str(iwavy)+" # IWAVY\n")
fid.write("5.5E-2 4 1 # ERVFILT, NFILT, IFILT\n")
fid.write("1 1 # IPA, NTH\n")
#fid.write(str(50.0*periodp_scl)+" 0.04 # TIMEP, TCP\n")
fid.write("2.0 0.04 # TIMEP, TCP\n")
fid.write("1 # TIMEW\n")
fid.write("0.0012285 # RDGL\n")
fid.write(str(20*DT)+" # TIMETURB\n")
fid.write(str(AKA)+" "+str(NSWAVEX)+" # AKA, NSWAVE\n")
fid.write(str(NSWAVEX)+" # NSWAVEX\n")
fid.write(str(USTAR)+" "+str(U10)+" "+str(USS)+" # USTAR, U10, USS\n")
fid.write(str(GAMMA)+" "+str(F)+" "+str(G)+" "+str(PHI)+" # gamma, F, G, PHI\n")
fid.write(str(IST)+" 3.3 # IST WEBER\n")
fid.close()

## fort.12
fid = open("fort.12.template", 'w')
fid.write(str(NXMOD)+" "+str(NYMOD)+" "+str(NPW)+" "\
  +str(NCPU)+" # NXMOD, NYMOD, NPW, NCPU\n")
fid.write("0 # ISTART\n")
fid.write(str(ISWELL)+" # ISWELL\n")
fid.write(str(NTIME)+" "+str(NOUTD)+" "+str(NOUTC)\
  +" # NTIME, NOUTD, NOUTC\n")
fid.write(str(NTP)+" "+str(DT_HOS0)+" # NTP, DT\n")
fid.write(str(PEX)+" "+str(PEY)+" # PEX, PEY\n")
fid.write(str(NSWAVEX)+" # NSWAVEX\n")
fid.write(str(AKA)+" "+str(NSWELLX)+" "+str(NSWELLY)\
  +" # AKA, NSWELLX, NSWELLY\n")
fid.write(str(FR2)+" # FR2\n")
fid.write(str(USTAR)+" "+str(U10)+" "+str(USS)\
  +" # USTAR, U10, USS\n")
fid.write(str(GAMMA)+" "+str(F)+" "+str(G)+" "\
  +str(PHI)+" # GAMMA, F, G, PHI\n")
fid.write(str(IST)+" "+str(BOND)+" "+str(ICAP)\
  +" # IST, BOND, ICAP\n")
fid.close()

## fort.13
fid = open("fort.13", 'w')
fid.write(str(NFHOS)+" "+str(NTHOS)+" # NFHOS, NTHOS\n")
fid.close()

## z_grid.inp

#NZ2 = NZ - NZ1
L2 = L12 - L1
L3 = ZL - L12

neven = 4

dz1 = np.zeros(NZ1)
dz2 = np.zeros(NZ2)
dz3 = np.zeros(NZ3+1)

q_1 = np.power(RAT1, 1.0/(NZ1 - neven - 1))
a1_1 = L1 / ((1.0 - np.power(q_1, NZ1-neven))/(1.0 - q_1) + 4.0)

for k in range(neven):
  dz1[k] = a1_1

for k in range(neven, NZ1):
  dz1[k] = a1_1 * np.power(q_1, k-neven)

print('dz1')
print(dz1)

q_3 = np.power(RAT3, 1.0/(NZ3-neven-1))
a1_3 = L3 / ((1.0 - np.power(q_3, NZ3-neven))/(1.0 - q_3) + 4.0*np.power(q_3, NZ3-neven-1))
for k in range(NZ3-neven):
  dz3[k] = a1_3*np.power(q_3,k) 

for k in range(NZ3-neven, NZ3):
  dz3[k] = dz3[NZ3-neven-1]


q_2 = np.power(RAT2, 1.0/(NZ2-1))
a1_2 = L2 / ((1.0 - np.power(q_2, NZ2))/(1.0 - q_2) )
for k in range(NZ2):
  dz2[k] = a1_2 * np.power(q_2, k) 
print('dz2')
print(dz2)

print('dz3')
print(dz3)

print('q_1='+str(q_1)+', a1_1='+str(a1_1))
print('q_2='+str(q_2)+', a1_2='+str(a1_2))
print('q_3='+str(q_3)+', a1_3='+str(a1_3))

dz12 = np.concatenate((dz1,dz2))
dz = np.concatenate((dz12,dz3))

dz1 = dz.copy()

dz[0] = dz1[0] / 2.
dz[NZ-2] = dz1[NZ-2] / 2.
dz[NZ-1] = dz1[NZ-1] / 2.

zz1 = np.zeros(NZ+1)
zz2 = np.zeros(NZ+1)
zz = np.zeros(NZ+1)
zz1[0] = 0.0
zz2[0] = 0.0

for k in range(1,NZ+1):
  zz1[k] = zz1[k-1] + dz1[k-1]
  zz2[k] = zz2[k-1] + dz[k-1]

for k in range(NZ+1):
  zz[k] = (zz2[k]-zz2[0])/(zz2[NZ-1]-zz2[0])*ZL

for k in range(NZ):
  dz[k] = zz[k+1] - zz[k]
#xx = np.zeros(NZ+1)
#plt.plot(xx, zz1, '+')
#plt.show()

result = np.column_stack((dz, zz, dz1, zz1, zz2))

np.savetxt("z_grid.dat", result, fmt='%.6f') 

fid = open("z_grid.inp",'w')
fid.write(str(NZ1)+" "+str(L1)+" "+str(RAT1)+"  # NZ1, L1, RAT1\n")
fid.write(str(NZ2)+" "+str(L2)+" "+str(RAT2)+"  # NZ2, L2, RAT2\n")
fid.write(str(NZ3)+" "+str(L3)+" "+str(RAT3)+"  # NZ3, L3, RAT3\n")
fid.close() 

print("grid around inner layer:")
for k in range(NZ1-3,NZ1+4):
  print(str(dz[k]))
print("grid around blade tip:")
for k in range(NZ1+NZ2-3,NZ1+NZ2+4):
  print(str(dz[k]))

# check bottom grid y+
print("yplus(k=1)="+str(utao_air_o)+"*"+str(dz[0]*l0)+"/"+str(nu_air)+"="+str(utao_air_o*dz[0]*l0/nu_air))

## control.dat
IR_scale = 1.0 / uinfty_o
if layout==3:
  IR_scale = 1.0/7.572136386

if layout==3 or layout==5 or layout==6 or layout==7:
  NumberOfNacelle = 3 * nturbine
  NumNacellePerLoc = 3
else:
  NumberOfNacelle = 1
  NumNacellePerLoc = 1

fid = open("control.dat",'w')
fid.write(str(rotor_model)+" # rotor_model\n")
fid.write(str(nturbine)+" "+str(nfoiltype)+" 3 # NumberOfTurbines, num_foiltype, num_blade\n")
fid.write("1.67 # loc_refvel\n")
fid.write("0 1 # FixTurbineAngvel, FixTipspeedRatio\n")
fid.write("1.0 # halfwidth_dfunc\n")
fid.write(str(nacelle_model)+" # nacelle_model\n")
fid.write("1 # rotate_nacelle\n")
fid.write(str(NumberOfNacelle)+" "+str(NumNacellePerLoc)+" # NumberOfNacelle, NumNacellePerLoc\n")
fid.write("1 1 # deltafunc_U, deltafunc_F\n")
fid.write(str(Shen1_AL)+" "+str(Shen1_AL_tipcorrection)+" "+str(Shen1_AL_tipcorrectionratio_Fa) \
  +"  # Shen1_AL, Shen1_AL_tipcorrection, Shen1_AL_tipcorrectionratio_Fa\n")
fid.write(str(InletRelaxation)+" "+str(IR_start)+" "+str(IR_end)+" "+str(IR_cx) \
  +" "+str(IR_cy)+" "+str(IR_cz)+" "+str(IR_NSmooth)+" "+str(IR_dataindex) \
  +" "+str(IR_scale) \
  +" # InletRelaxation IR_start IR_end IR_cx IR_cy IR_cz IR_NSmooth" \
  +" IR_dataindex IR_scale\n")
fid.write(str(c_bforce)+" # c_bforce\n")
fid.write(str(adjust_CT_ACL)+" "+str(adjust_CP_ACL)+" # adjust_CT_ACL, adjust_CP_ACL\n")
fid.write(str(adjust_CT_nacelle)+" # adjust_CT_nacelle\n")
fid.write(str(fsitype)+' # fsitype\n')
fid.close()

## Turbine.inp or turbine_dy_param.dat
if rotor_model==3 or rotor_model==5 or rotor_model==7:
  xij=np.zeros((nturbine, nturbine))
  yij=np.zeros((nturbine, nturbine))
  if layout==1 or layout==5 or layout==6 or layout==7:
    xij=np.zeros((nturbine_x, nturbine_y))
    yij=np.zeros((nturbine_x, nturbine_y))
    
    for i in range(nturbine_x):
      for j in range(nturbine_y):
        if layout!=7:
          xij[i,j]=(3.5+(nturbine_x/2.0-0.5)*7.0*np.sin(spacing_incl-inflowdir)-j*7.0*np.sin(spacing_incl-inflowdir)+i*7.0*np.cos(inflowdir))*d_scl
        elif layout==7:
          xij[i,j]=(3.5+(nturbine_x/2.0-0.5)*7.0*np.sin(spacing_incl-inflowdir)-j*7.0*np.sin(spacing_incl-inflowdir)+i*7.0*np.cos(inflowdir))*d_scl
        yij[i,j]=(-0.5*(nturbine_y-1.0)+j)*spacing_y_o/l0 + domain_y_scl/2.0
        print('i,j=('+str(i)+','+str(j)+'), (x,y)=('+str(xij[i,j])+','+str(yij[i,j])+')')

  if layout==2:
    xij[0,0]=0.6*domain_x_scl
    yij[0,0]=0.5*domain_y_scl

  if layout==3:
    xij[0,0]=2.5*d_scl
    yij[0,0]=0.5*domain_y_scl


  fid = open("Turbine.inp",'w')

  tmp="# nx_tb ny_tb nz_tb x_c y_c z_c Tipspeedratio angvel_fixed r_rotor pitch0 CT0 CP0\n"
  fid.write(tmp)

  for i in range(nturbine_x):
    for j in range(nturbine_y):
      tmp="1.0 0.0 0.0 "+str(xij[i,j])+" "+str(yij[i,j])+" "+str(hhub_scl)+" "+str(tsr_o)+" "+str(omega_scl)+" "+str(radius_scl)+" "+str(pitch_o)+" "+str(CT0_ACL)+" "+str(CP0_ACL)+"\n"
      fid.write(tmp)

  fid.close()

  fid = open("Turbine2.inp",'w')
  for i in range(nturbine_x):
    for j in range(nturbine_y):
      tmp = "0.0 0.0 "+str(omega_scl)+"\n"
      fid.write(tmp)

  fid.close()

  fid = open("FSI_Prescribe.inp", 'w')
  fsirefpos_scl = fsirefpos_o.copy()
  for i in range(nturbine_x):
    for j in range(nturbine_y):
      fsirefpos_scl[i,j,0] = xij[i,j] + 6.0/l0  # 6.0 is from rotor cener to tower center
      fsirefpos_scl[i,j,1] = yij[i,j]
        
  
  fsiamp_scl = fsiamp_o.copy()
  fsiamp_scl[0:3] = fsiamp_o[0:3]/l0
  fsiomega_scl = fsiomega_o.copy()
  fsiomega_scl = fsiomega_o * l0/ u0
  print(fsiomega_scl)
  print(fsiphase0)
  
  for i in range(nturbine_x):
    for j in range(nturbine_y):
      for k in range(6):
        fid.write(str(fsidof[k])+' '+str(fsirefpos_scl[i,j,k])+' '+str(fsiamp_scl[k])\
          +' '+str(fsiomega_scl[k])+' '+str(fsiphase0[k])+'\n')
  fid.close()

  fid = open("FSI2_rotor.inp", 'w')
  for i in range(nturbine_x):
    for j in range(nturbine_y):
      fid.write(str(xij[i,j])+' '+str(yij[i,j])+' '+str(hhub_scl)+' '\
        +'1.0 0.0 0.0\n')
  fid.close()

  fid = open("FSI2_nacelle.inp", 'w')
  for i in range(nturbine_x):
    for j in range(nturbine_y):
      for k in range(NumNacellePerLoc):
        fid.write(str(xij[i,j])+' '+str(yij[i,j])+' '+str(hhub_scl)+' '\
          +'1.0 0.0 0.0\n')
  fid.close()


elif rotor_model==1:
  fid = open("turbine_dy_param.inp",'w')
  fid.write(str(nturbine_x)+' '+str(nturbine_y)+' # nxwt, nywt\n')
  fid.write(str(radius_scl)+' '+str(hhub_scl)+' # rdisc, hdisc\n')
  fid.write(str(c_bforce)+" # c_bforce\n")
  fid.close()

## Nacelle.inp
if nacelle_model != 0:
  fid = open("Nacelle.inp", 'w')

  tmp = "# nx_tb ny_tb nz_tb x_c y_c z_c angvel_axis rotate_alongaxis frictionfactor dh xnacelle_upstreamend ynacelle_upstreamend znacelle_upstreamend CT0 area_CT ratio_CT0\n"
  fid.write(tmp)

  angvel_nac = np.zeros(NumNacellePerLoc)
  rotate_nac = np.zeros(NumNacellePerLoc, 'int')
  area_CT = np.zeros(NumNacellePerLoc)
  CT0_nac = np.zeros(NumNacellePerLoc)

  for k in range(NumNacellePerLoc):
    angvel_nac[k] = 0.0
    rotate_nac[k] = 0
    CT0_nac[k] = CT0_nacelle

  if NumNacellePerLoc==3:
    angvel_nac[1] = omega_scl
    rotate_nac[1] = 1
    CT0_nac[0] = 1.0
    CT0_nac[1] = 1.0
    CT0_nac[2] = 0.5


  if layout==5 or layout==6 or layout==7:
    area_CT[0] = 28.3309
    area_CT[1] = 28.3309
    area_CT[2] = 144.2
    ratio_CT0[0] = 0.00373
    ratio_CT0[1] = 0.01170
    ratio_CT0[2] = 0.00553

  for i in range(nturbine_x):
    for j in range(nturbine_y):
      for k in range(NumNacellePerLoc):
        tmp = "1.0 0.0 0.0 "+str(xij[i,j])+" "+str(yij[i,j])+" "+str(hhub_scl)+" "+str(angvel_nac[k])+" "\
          +str(rotate_nac[k])+" 0.0 "+str(np.power(dx_scl*dy_scl*dz_scl,1.0/3.0))+" -0.3 0.0 0.0 "\
          +str(CT0_nac[k])+" "+str(area_CT[k])+" "+str(ratio_CT0[k])+"\n"
        fid.write(tmp)

  fid.close()

## FOIL00
if layout==1 or layout==5 or layout==6 or layout==7:
#radius = 0.448799
  chord_orig_file = "chord.orig.dat"
  twist_orig_file = "twist.orig.dat"

  chord_orig = np.genfromtxt(chord_orig_file, comments="#")
  twist_orig = np.genfromtxt(twist_orig_file, comments="#")
# np.set_printoptions(precision=5)
#print(chord_orig)

  chord_npshape = chord_orig.shape
  if chord_orig[chord_npshape[0]-1,0]!=1.0 :
    temp = np.array([1.0, 0.0])
    chord_1 = np.row_stack((chord_orig, temp))

  twist_npshape = twist_orig.shape
  if twist_orig[twist_npshape[0]-1,0]!=1.0 :
    temp = np.array([1.0, 0.0])
    twist_1 = np.row_stack((twist_orig, temp))

  chord_npshape1 = chord_1.shape
  twist_npshape1 = twist_1.shape
  newsize = max(chord_npshape1[0],twist_npshape1[0])
  newstart = max(chord_1[0,0], twist_1[0,0])

  if chord_npshape1[0] >= twist_npshape1[0] :
    r1 = chord_1[:,0]
  else :
    r1 = twist_1[:,0] 

#print("before")
#print(r1)
  r1[0] = newstart
#print("middle:")
#print(r1)
  chord_2 = np.interp(r1, chord_1[:,0], chord_1[:,1])
  twist_2 = np.interp(r1, twist_1[:,0], twist_1[:,1])

  result_1 = np.column_stack((r1,twist_2,chord_2))
#print("result_1:")
#print(result_1)

  result_1[:,0] = result_1[:,0] * radius_scl
  result_1[:,2] = result_1[:,2] / 40.0 * radius_scl
#print("Rescaled:")
#print(result_1)

  h1 = "# Turbine type: recompiled Vestas V80 for Horns Rev 1\n# Airfoil type:\n"+str(newsize)
  np.savetxt("FOIL00",result_1,fmt='%.6f',header=h1,comments='')

if layout==2:
  nfoiltype = 5
  
  fname_all = 'FOIL_all'
  fall = open(fname_all,'w')
  for ifl in range(nfoiltype):
    fname_old = '{}{:02d}{}'.format('FOIL',ifl,'.orig')
    fname_new = '{}{:02d}'.format('FOIL',ifl)
    #print(fname_temp)
    fold = open(fname_old,'r')
    fnew = open(fname_new,'w')
    content = fold.readlines()
    #print(content)
    nline=int(str.strip(content[2]))
    #print(nline)
    for iline in range(nline):
      values = content[3+iline].split()
      #print(values)
      values[0] = str(float(values[0])/96.0/l0)
      values[2] = str(float(values[2])/96.0/l0)
      content[3+iline]=values[0]+" "+values[1]+" "+values[2]+"\n"
    fnew.writelines(content)
    fall.writelines(content)
    fold.close()
    fnew.close()

  fall.close()


## part2: acldata000
#if (layout==2 or layout==3 or layout==5):
#  newsize = np.floor((radius_scl-rmin_o/l0)/grid_res)+1
#  newsize = int(np.maximum(newsize,40))
#  r1 = np.linspace(rmin_o*2.0/d_o,1.0,newsize)

newsize = np.floor((radius_scl-rmin_o/l0)/grid_res)+1
newsize = int(np.maximum(newsize,40))
r1 = np.linspace(rmin_o*2.0/d_o,1.0,newsize)


fid = open("acldata000",'w')

nseg_acl = np.floor(radius_scl*(r1[newsize-1]-r1[0])/grid_res)
print("r1="+str(radius_scl*r1[0])+"~"+str(radius_scl*r1[newsize-1])+", grid_res="+str(grid_res)+", nseg_acl="+str(nseg_acl))

if nseg_acl <= newsize: 
  print("nseg_acl <= newsize: don't need refine the acldata000")
else :
  print("nseg_acl > newsize: need to refine acldata000 grid")

for ibi in range(3):
  tmp = str(newsize)+"\n"
  fid.write(tmp)
  # take care the sequence of blade index
  # make sure color (blade index) in acldata and acsdata are the same
  theta =  ibi * np.pi * 2.0/3.0
  for x in r1:
    xx1 = 0.0
    yy1 = x * radius_scl
    zz1 = 0.0
    xx2 = xx1
    yy2 = yy1 * np.cos(theta) + zz1 * np.sin(theta)
    zz2 = zz1 * np.cos(theta) - yy1 * np.sin(theta)
    tmp = str(xx2) + " " + str(yy2) + " " + str(zz2) + "\n"
    fid.write(tmp)

fid.close()


