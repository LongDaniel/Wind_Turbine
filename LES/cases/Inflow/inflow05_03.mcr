#!MC 1400
# Created by Tecplot 360 build 14.0.2.33360
#$!VarSet |MFBD| = 'C:\Users\lixx5444\Documents\plyu\PA\inflow05\'
$!VarSet |MFBD| = '/home/plyu/hosmnt/projects/PorteAgel_GWS_EP_6030/inflow05_2/'
#$!VarSet |MFBD| = '/home/pin/LES_HOS_Turbine/new_source_codes/note_LES/runjob/'

$!VarSet |NumLoop|=200
$!Loop |NumLoop|
$!VarSet |NewNo|=(|Loop|+200)
$!VarSet |Num|=(|NewNo|*250)

$!VarSet |toverperiod|=(|Num|*0.001)

$!READDATASET  '"-F" "1" "|MFBD|DAT_|Num%010d|.h5" "-D" "5" "/pp" "/u" "/v" "/w" "/z" "-I" "-K" "1" "1" "1"'
  DATASETREADER = 'HDF5 Loader'
$!ALTERDATA 
  EQUATION = '{x1}=({Z}-0.5)*6.2831853071795862/6.82954924693/128'
$!ALTERDATA 
  EQUATION = '{y1}=({Y}-0.5)*6.2831853071795862/8.72664625997/64'
$!ALTERDATA 
  EQUATION = '{Z}={/z}'
$!ALTERDATA 
  EQUATION = '{X}={x1}'
$!ALTERDATA 
  EQUATION = '{Y}={y1}'
$!ALTERDATA 
  EQUATION = '{Q} = 0'
$!GLOBALTHREEDVECTOR UVAR = 5
$!GLOBALTHREEDVECTOR VVAR = 6
$!GLOBALTHREEDVECTOR WVAR = 7
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
$!DELETEVARS  [12-35]
$!GLOBALCONTOUR 1  VAR = 4
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 1  LEGEND{SHOW = NO}

$!GLOBALCONTOUR 2  VAR = 5
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 2
  APPROXNUMVALUES = 15
  
$!CONTOURLEVELS NEW
  CONTOURGROUP = 2
  RAWDATA
11
0.00
0.10
0.20
0.30
0.40
0.50
0.60
0.70
0.80
0.90
1.00
$!GLOBALCONTOUR 2  COLORMAPFILTER{CONTINUOUSCOLOR{CMIN=0 CMAX = 1}}  

$!GLOBALCONTOUR 3  VAR = 6
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 3
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 4  VAR = 7
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 4
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 5  VAR = 11
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 5
  APPROXNUMVALUES = 15
 
$!GLOBALCONTOUR 2  COLORMAPFILTER{COLORMAPDISTRIBUTION = CONTINUOUS}
$!GLOBALCOLORMAP 1  CONTOURCOLORMAP = DIVBUYLRD

$!GLOBALCONTOUR 2  LEGEND{ISVERTICAL = NO}
$!GLOBALCONTOUR 2  LEGEND{ANCHORALIGNMENT = TOPCENTER}
$!GLOBALCONTOUR 2  LEGEND{XYPOS{X = 72}}
$!GLOBALCONTOUR 2  LEGEND{XYPOS{Y = 8}}
$!GLOBALCONTOUR 2  LEGEND{SHOW = YES}

$!ATTACHTEXT 
  ANCHORPOS
    {
    X = 5.0
    Y = 7.598784194528873
    }
  TEXT = 't=|toverperiod%.2f|'
  
$!THREEDAXIS FRAMEAXIS{XYPOS{X = 93}}
$!THREEDAXIS FRAMEAXIS{XYPOS{Y = 87.5}}

$!REDRAW 

$!VIEW FIT
$!SLICELAYERS SHOW = YES

$!SLICEATTRIBUTES 2  SHOWGROUP = YES
$!SLICEATTRIBUTES 2  SLICESURFACE = JPLANES

$!SLICEATTRIBUTES 1  SLICESURFACE = KPLANES

$!ISOSURFACELAYERS SHOW = YES
$!ISOSURFACEATTRIBUTES 1  DEFINITIONCONTOURGROUP = 5

$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 20

$!SLICEATTRIBUTES 1  CONTOUR{FLOODCOLORING = GROUP2}
$!SLICEATTRIBUTES 2  CONTOUR{FLOODCOLORING = GROUP2}
$!REDRAWALL 
$!ISOSURFACEATTRIBUTES 1  CONTOUR{FLOODCOLORING = GROUP2}

$!ISOSURFACEATTRIBUTES 1  EFFECTS{USETRANSLUCENCY = YES}

$!ROTATE3DVIEW THETA
  ANGLE = -30
  ROTATEORIGINLOCATION = DEFINEDORIGIN

$!VIEW FIT
$!REDRAWALL

$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 1200
$!EXPORTSETUP EXPORTFNAME = '|MFBD|export/export_style_3_ti_|NewNo%05d|.png'
$!EXPORT 
  EXPORTREGION = ALLFRAMES

$!EndLoop
$!RemoveVar |MFBD|