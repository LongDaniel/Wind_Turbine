#!MC 1400
# Created by Tecplot 360 build 14.0.2.33360
$!VarSet |MFBD| = '/home/plyu/hosmnt/projects/acs11'
$!READDATASET  '"-F" "1" "DAT_0000012000.h5" "-D" "7" "/pp" "/u" "/v" "/w" "/x" "/y" "/z" "-K" "1" "1" "1"'
  DATASETREADER = 'HDF5 Loader'
$!THREEDAXIS XDETAIL{VARNUM = 5}
$!THREEDAXIS YDETAIL{VARNUM = 6}
$!THREEDAXIS ZDETAIL{VARNUM = 7}
$!REDRAW 
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'SetFluidProperties Incompressible=\'T\' Density=1 SpecificHeat=2.5 UseSpecificHeatVar=\'F\' SpecificHeatVar=1 GasConstant=1 UseGasConstantVar=\'F\' GasConstantVar=1 Gamma=1.4 UseGammaVar=\'F\' GammaVar=1 Viscosity=1e-005 UseViscosityVar=\'F\' ViscosityVar=1 Conductivity=1 UseConductivityVar=\'F\' ConductivityVar=1'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'SetFieldVariables ConvectionVarsAreMomentum=\'F\' UVar=2 VVar=3 WVar=4 ID1=\'Pressure\' Variable1=1 ID2=\'NotUsed\' Variable2=0'
$!EXTENDEDCOMMAND 
  COMMANDPROCESSORID = 'CFDAnalyzer4'
  COMMAND = 'Calculate Function=\'SWIRL\' Normalization=\'None\' ValueLocation=\'CellCentered\' CalculateOnDemand=\'T\' UseMorePointsForFEGradientCalculations=\'F\''
$!GLOBALCONTOUR 1  VAR = 1
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!GLOBALCONTOUR 1  VAR = 8
$!CONTOURLEVELS RESETTONICE
  CONTOURGROUP = 1
  APPROXNUMVALUES = 15
$!CONTOURLEVELS NEW
  CONTOURGROUP = 1
  RAWDATA
11
-1
-0.9
-0.8
-0.7
-0.6
-0.5
-0.4
-0.3
-0.2
-0.1
0
$!SLICEATTRIBUTES 1  PRIMARYPOSITION{X = 0.78500000000000003}
$!SLICELAYERS SHOW = YES
$!REDRAW 
$!SLICEATTRIBUTES 1  SHOWGROUP = NO
$!SLICEATTRIBUTES 3  SHOWGROUP = YES
$!SLICEATTRIBUTES 3  PRIMARYPOSITION{Z = 0.78500000000000003}
$!REDRAW 
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 5
$!SLICELAYERS SHOW = NO
$!ISOSURFACELAYERS SHOW = YES
$!REDRAW 
$!ROTATE3DVIEW THETA
  ANGLE = -30
  ROTATEORIGINLOCATION = DEFINEDORIGIN
$!REDRAW 
$!ISOSURFACEATTRIBUTES 1  ISOVALUE1 = 3
$!REDRAW 
$!EXPORTSETUP EXPORTFORMAT = PNG
$!EXPORTSETUP IMAGEWIDTH = 1200
$!EXPORTSETUP EXPORTFNAME = 'swirl_iso_5_ti_12000.png'
$!EXPORT 
  EXPORTREGION = CURRENTFRAME
$!SAVELAYOUT  "swirl_iso_5.lay"
  USERELATIVEPATHS = YES
$!Quit
$!RemoveVar |MFBD|
