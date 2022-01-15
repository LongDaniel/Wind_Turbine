#!/usr/bin/python

'''
This is used for mapping the restart files of HOS-LES from an old case to a new case. 

Python Classes Objetcs: 
https://www.tutorialspoint.com/python/python_classes_objects.htm
http://thepythonguru.com/python-inheritance-and-polymorphism/
'''


import sys
import numpy as np
import scipy as sp
import scipy.interpolate as inp
import h5py

#class Grid(object):   # for python 2
class Grid():  # for python 3
  dim = 3
  g_n = np.ones((3), dtype='int') # (nz, ny, nx)
  g_v = np.ones((3,2)) # range of each dimension, first index is dimension, second is min and max. 
  g_L = np.ones((3)) # length of each dimension
  g_isUniform = np.ones((3), dtype='int')
  g_d = np.ones((3)) # dx of each dimension
  idbg = 1
  
  def __init__(self, ddim, nn, vv, isU, iidbg):   
    self.dim = ddim
    self.g_n[0:dim] = nn
    self.g_v[0:dim, :] = vv
    self.g_isUniform[0:dim] = isU
    self.idbg = iidbg
    for i in range(self.dim):
      self.g_L[i]= self.g_v[i,1]  - self.g_v[i,0]
      if self.g_isUniform[i]==1:
        self.g_d[i] = self.g_L[i] / self.g_n[i]
      if self.idbg==1:
        print("Grid in {:d}th direction: n={:d}, range=[{:e5},{:e5}], isUniform={:d}".format(i, self.g_n[i], self.g_v[i,0], self.g_v[i,1], self.g_isUniform[i]))
      
  
class HosLesGrid(Grid):
  idbg = 1
  dir_case = ""
  hlg_fn = ""
  zz = []
  zw = []
  x1d = [] # x value index from 1 to nx
  y1d = [] 
  x1d2 = [] # x value index from 0 to nx, so that it will not exceed input range during the interpolation
  y1d2 = []

  x1d_h = []
  y1d_h = []
  x1d2_h = [] # for hos grid
  y1d2_h = []
  
  coord2d_xy = []
  coord3d = []
  coord3d_w = []
  
  
  def __init__(self, ddim, nn, vv, isU, iidbg):
    #super(HosLesGrid, self).__init__(ddim, nn, vv, isU, iidbg) # for python 2
    super().__init__(ddim, nn, vv, isU, iidbg) # for python 3
    self.idbg = iidbg 
    
  def __init__(self, dircase):
    self.dir_case = dircase
    print("Initializing HosLesGrid from dir: {:s}".format(self.dir_case))
    self.hlg_fn = dircase + "/grid.h5"
    #self.set_z_from_h5(self.hlg_fn)
  
  def set_z(self, zzin, zwin):
    self.zz = zzin
    self.zw = zwin
    
  def set_z_from_h5(self, fname):
    print("Reading z grid from H5 file: "+fname)
    f = h5py.File(fname, "r")
    self.zz = np.array(f["zz"])
    self.zw = np.array(f["zw"])
    f.close()   
    if self.idbg>0:
      print("zz readed. [{:f}, {:f}] (len={:d})".format(self.zz[0], self.zz[-1], len(self.zz)))
      print("zw readed. [{:f}, {:f}] (len={:d})".format(self.zw[0], self.zw[-1], len(self.zw)))
    
class HosLesField(HosLesGrid):
  idbg = 1
  u = []
  v = []
  w = []
  p = []
  eta = []
  hh = []
  tmp3d = []
  tmp3d2 = []
  tmp2d = []
  tmp2d2 = []
  hlf_fn_les_h5 = ""
  hlf_fn_les_aux = ""
  hlf_fn_les_param = ""
  hlf_fn_hos_h5 = ""
  hlf_fn_hos_param = ""
  ihos = 0
  hg = [] # hos_grid
   
  def __init__(self, dircase, pex, pey, hbar, ihos = 0):
    #super(HosLesField, self).__init__(dircase) # for python 2
    super().__init__(dircase) # for python 3
    super().set_z_from_h5(self.hlg_fn)
    print("Initializing HosLesField from dir: {:s}".format(self.dir_case))
    self.g_L = np.array((hbar, np.pi * 2.0 / pey, np.pi * 2.0 / pex))
    self.g_isUniform = np.array((0, 1, 1), dtype='int')
    self.hlf_fn_les_h5 = dircase + "/restart.h5"  
    self.hlf_fn_les_aux = dircase + "/restart_aux.h5" 
    self.ihos = ihos
    self.hlf_fn_hos_h5 = dircase + "/restart_hos.h5"  
    self.read_hlf_var(fi=1, mode=0)
    if ihos == 1:
      self.hg = HosLesGrid(dircase)
      self.read_hlf_var(fi=3, mode=0)
    
  def read_hlf_var(self, fi=1, mode=1, key="u"):
    if fi==1:
      fn = self.hlf_fn_les_h5
    elif fi==2:
      fn = self.hlf_fn_les_aux
    elif fi==3:
      fn = self.hlf_fn_hos_h5
    
    print("Reading HosLesField data, filename={:s}".format(fn))
    f = h5py.File(fn, "r")    
    if mode==0 and fi==1: # test mode, to get the shape of les data
      self.tmp3d = f["u"]
      self.g_n = np.array(self.tmp3d.shape) # [nz,ny,nx]
      self.g_d = np.array((-1.0, self.g_L[1]/self.g_n[1], self.g_L[2]/self.g_n[2]))
      
      self.x1d = self.g_d[2] * np.linspace(1, self.g_n[2], num=self.g_n[2])
      self.y1d = self.g_d[1] * np.linspace(1, self.g_n[1], num=self.g_n[1])
      self.x1d2 = self.g_d[2] * np.linspace(0, self.g_n[2], num=self.g_n[2]+1)
      self.y1d2 = self.g_d[1] * np.linspace(0, self.g_n[1], num=self.g_n[1]+1)
      
      print("Readed HosLesField les data in test mode, 3d_shape={:s}\n".format(str(self.g_n)))
    
    elif mode==0 and fi==3 and self.ihos==1: # test mode, to get the shape of les data
      self.tmp2d = f["eta_hos"]
      self.hg.g_n = np.array([-1, self.tmp2d.shape[0], self.tmp2d.shape[1]]) # [-1,ny,nx]
      self.hg.g_d = np.array([-2.0, self.g_L[1]/self.hg.g_n[1], self.g_L[2]/self.hg.g_n[2]])
      self.x1d_h = self.hg.g_d[2] * np.linspace(1, self.hg.g_n[2], num=self.hg.g_n[2])
      self.y1d_h = self.hg.g_d[1] * np.linspace(1, self.hg.g_n[1], num=self.hg.g_n[1])
      self.x1d2_h = self.hg.g_d[2] * np.linspace(0, self.hg.g_n[2], num=self.hg.g_n[2]+1)
      self.y1d2_h = self.hg.g_d[1] * np.linspace(0, self.hg.g_n[1], num=self.hg.g_n[1]+1)
    
      print("Readed HosLesField hos data in test mode, 2d_shape={:s}\n".format(str(self.hg.g_n[1:3])))
      
      #for i in range(self.g_n[2]):
      # for j in range(self.g_n[1]):
      #   self.coord2d_xy[j, i] = 
      
      #self.coord3d = np.zeros(self.g_n)
      #self.coord3d_w = np.zeros(self.g_n)
      #for k in range(self.g_n[0]):
      # self.coord3d[]
      
    elif mode==1 or mode==2: # read only one 3D variable with the key (default mode)
      self.tmp3d = np.array(f[key])
      n = self.g_n
      self.tmp3d2 = np.empty((n[0], n[1]+1, n[2]+1))
      self.tmp3d2[:,0,1:] = self.tmp3d[:,n[1]-1,:]
      self.tmp3d2[:,1:,0] = self.tmp3d[:,:,n[2]-1]
      self.tmp3d2[:,1:,1:] = self.tmp3d
      self.tmp3d2[:,0,0] = self.tmp3d2[0,n[1]-1,0]      
      print("Readed HosLesField les data '{:s}', 3d_shape={:s}".format(key, str(self.tmp3d.shape)))
    elif mode==3: # read an 2D variable in xy plane with the key
      self.tmp2d = np.array(f[key])
      if fi == 3:
        n = self.hg.g_n
      else:
        n = self.g_n
      self.tmp2d2 = np.empty((n[1]+1, n[2]+1))
      self.tmp2d2[0,1:] = self.tmp2d[n[1]-1,:]
      self.tmp2d2[1:,0] = self.tmp2d[:,n[2]-1]
      self.tmp2d2[1:,1:] = self.tmp2d
      self.tmp2d2[0,0] = self.tmp2d2[n[1]-1,0]
      if fi == 3:
        print("Readed HosLesField hos data '{:s}', 2d_shape={:s}".format(key, str(self.tmp2d.shape))) 
      else:
        print("Readed HosLesField les data '{:s}', 2d_shape={:s}".format(key, str(self.tmp2d.shape))) 

    f.close()
    
  #def __del__(self):
  # f.close()
  
def interpOneVariable(old, fnew, keyin, pts, fiin, modein):
  print("Begin interpolation for key='{:s}', newsize={:s}".format(keyin, str(fnew[keyin].shape)))
  if modein==1: # 3D variables at zz
    old.read_hlf_var(fi=fiin, key=keyin, mode=1)
    print("grid shape=({:d}, {:d}, {:d}), value shape={:s}".format(len(old.zz), len(old.y1d2), len(old.x1d2), str(old.tmp3d2.shape)))
    print("The source domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}], z:[{4:f}, {5:f}]".format(\
      np.min(old.x1d2), np.max(old.x1d2), np.min(old.y1d2), np.max(old.y1d2), \
      np.min(old.zz), np.max(old.zz) ))
    interp_func = inp.RegularGridInterpolator((old.zz, old.y1d2, old.x1d2), old.tmp3d2) 
    v = fnew[keyin]
    v[:,:,:] = interp_func(pts)
  elif modein==2: # 3D variables at zw
    old.read_hlf_var(fi=fiin, key=keyin, mode=modein)
    print("grid shape=({:d}, {:d}, {:d}), value shape={:s}".format(len(old.zw), len(old.y1d2), len(old.x1d2), str(old.tmp3d2.shape)))
    print("The source domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}], z:[{4:f}, {5:f}]".format(\
      np.min(old.x1d2), np.max(old.x1d2), np.min(old.y1d2), np.max(old.y1d2), \
      np.min(old.zw), np.max(old.zw) ))
    interp_func = inp.RegularGridInterpolator((old.zw, old.y1d2, old.x1d2), old.tmp3d2) 
    v = fnew[keyin]
    # the top grid zw might be out of range of old.zw, so ignore it at the first step.
    ss = pts.shape
    v[:(ss[0]-1),:,:] = interp_func(pts[:(ss[0]-1),:,:])
    v[ss[0]-1,:,:] = 2.0 * v[ss[0]-2,:,:] - v[ss[0]-3,:,:]
  elif modein==3: # 2D variable at xy plane
    old.read_hlf_var(fi=fiin, key=keyin, mode=modein)
    if fiin == 3:
      print("HOS grid shape=({:d}, {:d}), value shape={:s}".format(len(old.y1d2_h), len(old.x1d2_h), str(old.tmp2d2.shape)))
      print("The source domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}]".format(\
        np.min(old.x1d2_h), np.max(old.x1d2_h), np.min(old.y1d2_h), np.max(old.y1d2_h) )) 
      interp_func = inp.RegularGridInterpolator((old.y1d2_h, old.x1d2_h), old.tmp2d2) 
    else:
      print("LES grid shape=({:d}, {:d}), value shape={:s}".format(len(old.y1d2), len(old.x1d2), str(old.tmp2d2.shape)))
      print("The source domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}]".format(\
        np.min(old.x1d2), np.max(old.x1d2), np.min(old.y1d2), np.max(old.y1d2) )) 
      interp_func = inp.RegularGridInterpolator((old.y1d2, old.x1d2), old.tmp2d2) 
    v = fnew[keyin]
    v[:,:] = interp_func(pts)

  if modein == 3:
    print("The target domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}]".format(\
      np.min(pts[:,:,1]), np.max(pts[:,:,1]), np.min(pts[:,:,0]), np.max(pts[:,:,0])))
  else:
    print("The target domain is x:[{0:f}, {1:f}], y:[{2:f}, {3:f}], z:[{4:f}, {5:f}]".format(\
      np.min(pts[:,:,:,2]), np.max(pts[:,:,:,2]), np.min(pts[:,:,:,1]), np.max(pts[:,:,:,1]), \
      np.min(pts[:,:,:,0]), np.max(pts[:,:,:,0])))
  print("End interpolation for key='{:s}'\n".format(keyin))

def interpHosLesField(old, new):
  # Based on: https://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.interpolate.RegularGridInterpolator.html

  pts = np.empty((new.g_n[0], new.g_n[1], new.g_n[2], 3)) # g_n[0]*g_n[1]*g_n[2] size of 3-element vector
  pts_w = np.empty((new.g_n[0], new.g_n[1], new.g_n[2], 3))
  pts2d = np.empty((new.g_n[1], new.g_n[2], 2))
  if old.ihos==1 and new.ihos==1:
    pts2d_h = np.empty((new.hg.g_n[1], new.hg.g_n[2], 2)) # hos grid

  for i in range(new.g_n[2]):
    for j in range(new.g_n[1]):
      pts2d[j,i,:] = np.array(new.y1d[j], new.x1d[i])
      for k in range(new.g_n[0]):
        pt = np.array((new.zz[k], new.y1d[j], new.x1d[i]))
        pt_w = np.array((new.zw[k], new.y1d[j], new.x1d[i]))
        # print("tmp3d[{:d},{:d},{:d}]={:s}".format(k, j, i, str(interp_func(pt))))
        pts[k,j,i, :] = pt
        pts_w[k,j,i,:] = pt_w
  
  if old.ihos==1 and new.ihos==1:
    for i in range(new.hg.g_n[2]):
      for j in range(new.hg.g_n[1]):
        if old.ihos==1 and new.ihos==1:
          pts2d_h[j,i,:] = np.array(new.y1d_h[j], new.x1d_h[i])
  
  f = h5py.File(new.hlf_fn_les_h5, 'r+')
  interpOneVariable(old, f, "eta", pts2d, 1, 3)
  interpOneVariable(old, f, "hh", pts2d, 1, 3)
  interpOneVariable(old, f, "u", pts, 1, 1)
  interpOneVariable(old, f, "v", pts, 1, 1)
  interpOneVariable(old, f, "w", pts_w, 1, 2)
  interpOneVariable(old, f, "pp", pts, 1, 1)
  f.close()
  
  f = h5py.File(new.hlf_fn_les_aux, 'r+')
  interpOneVariable(old, f, "hu", pts, 2, 1)
  interpOneVariable(old, f, "hv", pts, 2, 1)
  interpOneVariable(old, f, "hw", pts, 2, 2)
  interpOneVariable(old, f, "ht", pts2d, 2, 3)
  interpOneVariable(old, f, "u1", pts2d, 2, 3)
  interpOneVariable(old, f, "v1", pts2d, 2, 3)
  f.close()
  
  if old.ihos==1 and new.ihos==1:
    f = h5py.File(new.hlf_fn_hos_h5, 'r+')
    interpOneVariable(old, f, "eta_hos", pts2d_h, 3, 3)
    interpOneVariable(old, f, "vps_hos", pts2d_h, 3, 3)
    interpOneVariable(old, f, "pa_hos", pts2d_h, 3, 3)
  
if len(sys.argv)!=3:
  print("Error: number of args incorrect.")
  print("Use it in this way: \"python pyMapField.py old_case_name new_case_name\"")
  exit(1)

print("numpy version="+np.__version__)
print("scipy version="+sp.__version__+" (should > 0.14)")

old_name = sys.argv[1]
new_name = sys.argv[2]
print("Reading args: old_name={:s}, new_name={:s}\n".format(old_name, new_name))

old_field = HosLesField(old_name, 1.45444104333, 8.72664625997, 0.46, ihos=0)
new_field = HosLesField(new_name, 1.45444104333, 8.72664625997, 0.46, ihos=0)
#old_field = HosLesField(old_name, 0.00280499344071, 0.00560998688141, 500.0, ihos=0)
#new_field = HosLesField(new_name, 0.00280499344071, 0.00560998688141, 500.0, ihos=0)
#old_field = HosLesField(old_name, 0.00280499344071, 0.00280499344071, 500.0, ihos=1)
#new_field = HosLesField(new_name, 0.00280499344071, 0.00280499344071, 500.0, ihos=1)

interpHosLesField(old_field, new_field)

print("Other files to be modified.")
print("fort.11: dt, pex, pey, nx, ny, nz")

cmd = "cp old/restart_param.dat new/"
print(cmd)
cmd = "cp old/restart_param_hos.dat new/"
print(cmd)
cmd = "cp old/restart_uref.dat new/"
print(cmd)
cmd = "cp old/Turbine2.inp new/"
print(cmd)
cmd = "cp old/FSI*.inp new/"
print(cmd)
