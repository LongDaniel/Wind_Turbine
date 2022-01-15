import numpy as np

radius = 0.448799
chord_orig_file = "chord.orig.dat"
twist_orig_file = "twist.orig.dat"

## part1: FOIL00

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
print("result_1:")
print(result_1)

result_1[:,0] = result_1[:,0] * radius
result_1[:,2] = result_1[:,2] / 40.0 * radius
print("Rescaled:")
print(result_1)

h1 = "# Turbine type: recompiled Vestas V80 for Horns Rev 1\n# Airfoil type:\n"+str(newsize)
np.savetxt("FOIL00",result_1,fmt='%.6f',header=h1,comments='')

## part2: acldata000

fid = open("acldata000",'w')

for ibi in range(3):
  tmp = str(newsize)+"\n"
  fid.write(tmp)
  # take care the sequence of blade index
  # make sure color (blade index) in acldata and acsdata are the same
  theta = - ibi * np.pi * 2.0/3.0
  for x in r1:
    xx1 = 0.0
    yy1 = x * radius
    zz1 = 0.0
    xx2 = xx1
    yy2 = yy1 * np.cos(theta) + zz1 * np.sin(theta)
    zz2 = zz1 * np.cos(theta) - yy1 * np.sin(theta)
    tmp = str(xx2) + " " + str(yy2) + " " + str(zz2) + "\n"
    fid.write(tmp)

fid.close()


### part3: Turbine.inp
#nturbine = 2 
#
#
#lref=7.0*radius/(2.0*np.pi)
#xij=np.zeros((nturbine, nturbine))
#yij=np.zeros((nturbine, nturbine))
#alp=7.0/180.0*np.pi
#for i in range(nturbine):
#  for j in range(nturbine):
#    xij[i,j]=(3.5+(nturbine/2.0-0.5)*7.0*np.sin(alp)-j*7.0*np.sin(alp)+i*7.0)*radius/lref
#    yij[i,j]=(3.5+7*j)*np.cos(alp)*radius/lref
#
#fid = open("Turbine.inp",'w')
#
#tmp="# nx_tb ny_tb nz_tb x_c y_c z_c Tipspeedratio angvel_fixed r_rotor pitch0\n"
#fid.write(tmp)
#
#for i in range(nturbine):
#  for j in range(nturbine):
#    tmp="1.0 0.0 0.0 "+str(xij[i,j])+" "+str(yij[i,j])+" 0.785398 8.45235 14.80394 "+str(radius)+" -1.02\n"
#    fid.write(tmp)
#
#fid.close()
