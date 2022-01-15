#!MC 1400
# Created by Tecplot 360 build 14.0.2.33360
$!VarSet |MFBD| = '/home/plyu/hosmnt/projects/acs11'
$!READDATASET  '"-F" "1" "DAT_0000016900.h5" "-D" "7" "/pp" "/u" "/v" "/w" "/x" "/y" "/z" "-K" "1" "1" "1"'
  DATASETREADER = 'HDF5 Loader'
$!THREEDAXIS XDETAIL{VARNUM = 5}
$!THREEDAXIS YDETAIL{VARNUM = 6}
$!THREEDAXIS ZDETAIL{VARNUM = 7}
$!ALTERDATA 
  EQUATION = '{Q} = 0'
$!GLOBALTHREEDVECTOR UVAR = 2
$!GLOBALTHREEDVECTOR VVAR = 3
$!GLOBALTHREEDVECTOR WVAR = 4
$!ALTERDATA 
  EQUATION = '{dudx} = ddx(u)'
$!ALTERDATA 
  EQUATION = '{dvdx} = ddx(v)'
$!ALTERDATA 
  EQUATION = '{dwdx} = ddx(w)'
$!ALTERDATA 
  EQUATION = '{dudy} = ddy(u)'
$!ALTERDATA 
  EQUATION = '{dvdy} = ddy(v)'
$!ALTERDATA 
  EQUATION = '{dwdy} = ddy(w)'
$!ALTERDATA 
  EQUATION = '{dudz} = ddz(u)'
$!ALTERDATA 
  EQUATION = '{dvdz} = ddz(v)'
$!ALTERDATA 
  EQUATION = '{dwdz} = ddz(w)'
$!ALTERDATA 
  EQUATION = '{s11} = {dudx}'
$!ALTERDATA 
  EQUATION = '{s12} = 0.5*({dudy}+{dvdx})'
$!ALTERDATA 
  EQUATION = '{s13} = 0.5*({dudz}+{dwdx})'
$!ALTERDATA 
  EQUATION = '{s22} = {dvdy}'
$!ALTERDATA 
  EQUATION = '{s23} = 0.5*({dvdz}+{dwdy})'
$!ALTERDATA 
  EQUATION = '{s33} = {dwdz}'
$!ALTERDATA 
  EQUATION = '{Omga12} = 0.5*({dudy}-{dvdx})'
$!ALTERDATA 
  EQUATION = '{Omga13} = 0.5*({dudz}-{dwdx})'
$!ALTERDATA 
  EQUATION = '{Omga23} = 0.5*({dvdz}-{dwdy})'
$!ALTERDATA 
  EQUATION = '{s2o2_11} = {s11}**2 + {s12}**2 + {s13}**2 - {Omga12}**2 - {Omga13}**2'
$!ALTERDATA 
  EQUATION = '{s2o2_12} = {s11}*{s12} + {s12}*{s22} + {s13}*{s23} - {Omga13}*{Omga23}'
$!ALTERDATA 
  EQUATION = '{s2o2_13} = {s11}*{s13} + {s12}*{s23} + {s13}*{s33} - {Omga12}*{Omga23}'
$!ALTERDATA 
  EQUATION = '{s2o2_22} = {s12}**2 + {s22}**2 + {s23}**2 - {Omga12}**2 - {Omga23}**2'
$!ALTERDATA 
  EQUATION = '{s2o2_23} = {s12}*{s13} + {s22}*{s23} + {s23}*{s33} - {Omga12}*{Omga13}'
$!ALTERDATA 
  EQUATION = '{s2o2_33} = {s13}**2 + {s23}**2 + {s33}**2 - {Omga13}**2 - {Omga23}**2'
$!ALTERDATA 
  EQUATION = '{Q} = 2*{Omga12}**2 + 2*{Omga13}**2 + 2*{Omga23}**2 - {S11}**2 - {S22}**2 - {S33}**2 - 2*{S12}**2 - 2*{S13}**2 - 2*{S23}**2'
$!DELETEVARS  [9-32]
$!GLOBALCONTOUR 1  VAR = 1
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 1  VAR = 8
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 12
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 15
$!ISOSURFACELAYERS SHOW = YES
$!REDRAW 
$!ROTATE3DVIEW THETA
  ANGLE = -120
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!REDRAW 
$!ROTATE3DVIEW THETA
  ANGLE = -30
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!REDRAW 
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 20
$!REDRAW 

$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 1200
$!EXPORTSETUP EXPORTFNAME = 'export_q_iso_20_ti_0000016900.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!Quit
$!RemoveVar |MFBD|
