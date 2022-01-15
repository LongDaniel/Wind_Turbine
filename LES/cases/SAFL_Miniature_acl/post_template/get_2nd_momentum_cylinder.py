
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# In[1]:


import numpy as np
import scipy as sp
import scipy.interpolate as inp

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import cmocean

import h5py


# In[2]:


def find_index(_z, _limits):
    _n = len(_z)
    _i_min = 0
    _i_max = _n - 1
    _limits2 = [0.0,0.0]
    if isinstance(_limits, float):
        _limits2[0] = _limits
        _limits2[1] = _limits
    else:
        _limits2 = _limits
        
    for i in range(_n):
        if _z[i]<_limits2[0] and i>_i_min :
            _i_min = i
        if _z[i]>_limits2[1] and i<_i_max :
            _i_max = i
            
    #print('zlimits='+str(_limits))
    #print('i_min='+str(_i_min)+', i_max='+str(_i_max))
    _i_min += 1
    
    if isinstance(_limits, float):
        return int(_i_min)
    else:
        return int(_i_min), int(_i_max)


# In[3]:


def get_var_cart(fn, mode = 0):
    if mode == 0:
        # Instantaneous field
        f = h5py.File(fn, 'r')
        u_cart = np.array(f['u'])
        v_cart = np.array(f['v'])
        w2_cart = np.array(f["w"])
        w_cart = np.empty(w2_cart.shape)
        
        ## interpolate w to u,v,p grid.
        for k in range(nn[2]):
            if k==0:
                w_cart[k,:,:] = w2_cart[k,:,:]
            else:
                w_cart[k,:,:] = 0.5*w2_cart[k,:,:]+0.5*w2_cart[k-1,:,:]

        f.close()
        
        d = {'u_cart': u_cart, 'v_cart': v_cart, 'w_cart': w_cart}
    elif mode == 1:
        # Mean field
        #fn = "mean_field_3d.h5"
        f = h5py.File(fn, 'r')
        #print(list(f.keys()))

        u_cart = np.array(f['u'])
        v_cart = np.array(f['v'])
        w2_cart = np.array(f["w"])
        w_cart = np.empty(w2_cart.shape)
        nut_cart = np.array(f['nut'])

        ## interpolate w to u,v,p grid.
        for k in range(nn[2]):
            if k==0:
                w_cart[k,:,:] = w2_cart[k,:,:]
            else:
                w_cart[k,:,:] = 0.5*w2_cart[k,:,:]+0.5*w2_cart[k-1,:,:]

        f.close()
        
        d = {'u_cart': u_cart, 'v_cart': v_cart, 'w_cart': w_cart, 'nut_cart': nut_cart}
    else:
        raise ValueError("Mode not implemented.")
    return d


# In[4]:


def cart_to_cyl_2d(_x, _y, center = [0,0]):
    _x2 = _x - center[0]
    _y2 = _y - center[1]
    _r = np.sqrt(_x2**2 + _y2**2)
    _theta = np.arctan2(_y2, _x2)
    return _r, _theta
def cyl_to_cart_2d(_r, _theta, center = [0,0]):
    _x = _r * np.cos(_theta) + center[0]
    _y = _r * np.sin(_theta) + center[1]
    return _x, _y


# In[5]:


def get_vel_cyl(_u, _v, _w, _x, _y, _z, _cylpts):
    # Interpolate the u, v, w to controlling grids in cylindrical coordinates
    interp_func = inp.RegularGridInterpolator((_x, _y, _z), _u)
    _u_cyl = interp_func(_cylpts)

    interp_func = inp.RegularGridInterpolator((_x, _y, _z), _v)
    _v_cyl = interp_func(_cylpts)

    interp_func = inp.RegularGridInterpolator((_x, _y, _z), _w)
    _w_cyl = interp_func(_cylpts)

    # Transform u, v, w in directions of cartisian axes to u, ur, utheta in directions of cylindrical axes
    _ur_cyl = np.empty(_v_cyl.shape)
    _uth_cyl = np.empty(_w_cyl.shape)
    
    for i in range(_cylpts.shape[0]):
        for j in range(_cylpts.shape[1]):
            for k in range(_cylpts.shape[2]):
                tmp_cos = np.cos(thth[k])
                tmp_sin = np.sin(thth[k])
                _ur_cyl[i,j,k] = _v_cyl[i,j,k] * tmp_cos + _w_cyl[i,j,k] * tmp_sin
                _uth_cyl[i,j,k] = _w_cyl[i,j,k] * tmp_cos - _v_cyl[i,j,k] * tmp_sin
    return _u_cyl, _ur_cyl, _uth_cyl


# In[25]:


def grad_1d(_f, _x):
    '''Return the first order derivative with second order accuracy.'''
    _n = len(_x)
    _dx = _x[1] - _x[0] 
    _dfdx = np.empty(_n)
    #for i in range(1, _n-1):
    #    _dfdx[i] = (_f[i+1] - _f[i-1]) / _dx / 2.0
    _dfdx[1:(_n-1)] = (_f[2:_n] - _f[0:(_n-2)]) / _dx / 2.0
    _dfdx[0] = (-1.5*_f[0]+2.0*_f[1]-0.5*_f[2]) / _dx
    _dfdx[_n-1] = (1.5*_f[_n-1]-2.0*_f[_n-2]+0.5*_f[_n-3]) / _dx
    return _dfdx

def grad_3d(_f3d, _x, axis = 0):
    _sz = _f3d.shape
    _n = len(_x)
    if _sz[axis] != _n:
        print('Error: size incorrect.')
    _df3d = np.empty(_sz)
    if axis == 0:
        for j in range(_sz[1]):
            for k in range(_sz[2]):
                _df3d[:,j,k] = grad_1d(_f3d[:,j,k], _x)
    elif axis == 1:
        for i in range(_sz[0]):
            for k in range(_sz[2]):
                _df3d[i,:,k] = grad_1d(_f3d[i,:,k], _x)
    elif axis == 2:
        for i in range(_sz[0]):
            for j in range(_sz[1]):
                _df3d[i,j,:] = grad_1d(_f3d[i,j,:], _x)
    return _df3d  

def test_grad_1d():
    test_x = np.linspace(0, 1, 100)
    test_y = test_x ** 2
    test_dydx = grad_1d(test_y, test_x)

    fig, ax = plt.subplots()
    ax.plot(test_x, test_dydx, '+-')
    plt.show()
    plt.close()
    return


# In[7]:


# timestep range
ti_s = 9030 
ti_e = 33000 + 1
ti_i = 30
ti_arr = np.arange(ti_s, ti_e, ti_i)
ti_n = len(ti_arr)

#foldername = "F:\\UMN\\remote_temp\\PA_validation\\"
foldername = "./"
fn_grid_h5 = 'grid.h5'

## domain length, grid number, grid cell size
ll = np.array([4.32, 0.72, 0.46])
nn = np.array([384, 96, 65])
dd = ll / nn

## grid
g = h5py.File(foldername+fn_grid_h5, "r")
zz = np.array(g["zz"])
zw = np.array(g["zw"])
g.close()
xx = dd[0] * np.arange(nn[0])
yy = dd[1] * np.arange(nn[1])
#print(len(zz))
#print(zz)

## turbine
loc_tbn = np.array([0.75, 0.36, 0.125])
r_tbn = 0.075


# In[8]:


center_cyl = loc_tbn[1:3]

rr_max = 1.5 * r_tbn
nn_r = 50
rr = rr_max / nn_r * np.arange(1, nn_r + 1)
#print('rr = \n'+str(rr))

nn_theta = 120
thth = np.pi * 2.0 / nn_theta * np.arange(1, nn_theta + 1)
#print('thth = \n'+str(thth))

xx2_min = loc_tbn[0] - 3.0 * r_tbn
xx2_max = loc_tbn[0] + 24.0 * r_tbn
#print([xx2_min, xx2_max])
imin_xx2, imax_xx2 = find_index (xx, [xx2_min, xx2_max])
nn_xx2 = imax_xx2 - imin_xx2
xx2 = xx[imin_xx2:imax_xx2]
#print(xx2)

ix_tb = find_index(xx2, loc_tbn[0])
#print('ix_tb in xx2 is '+str(ix_tb))

xx_D = (xx2 - loc_tbn[0]) / (2.0 * r_tbn)
rr_D = rr / (2.0 * r_tbn)


# In[9]:


#nn_cyl = np.array([nn_xx2, nn_r, nn_theta])
#u_cyl = np.empty(nn_cyl)
#v_cyl = np.empty(nn_cyl)
#w_cyl = np.empty(nn_cyl)
#nut_cyl = np.empty(nn_cyl)


# In[10]:


cartpts_cyl = np.empty([nn_xx2, nn_r, nn_theta, 3])
for i in range(nn_xx2):
    cartpts_cyl[i, :, :, 2] = xx2[i]
    for j in range(nn_r):
        for k in range(nn_theta):
            cartpts_cyl[i, j, k, 1], cartpts_cyl[i, j, k, 0] = cyl_to_cart_2d(rr[j], thth[k], center = center_cyl)


# In[11]:


#fn = 'mean_field_3d.h5'
#dict_vel = get_var_cart(fn, mode = 1)
#u_cart, v_cart, w_cart = dict_vel['u_cart'], dict_vel['v_cart'], dict_vel['w_cart']
#ux_cyl_m, ur_cyl_m, uth_cyl_m = get_vel_cyl(u_cart, v_cart, w_cart, zz*ll[2], yy, xx, cartpts_cyl)


# In[12]:


i_main = True 
if i_main:
    ux_cyl_m2 = np.zeros([nn_xx2, nn_r, nn_theta])
    ur_cyl_m2 = np.zeros([nn_xx2, nn_r, nn_theta])
    uth_cyl_m2 = np.zeros([nn_xx2, nn_r, nn_theta])
    uxpuxp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    urpurp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    uthputhp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    uxpurp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    uxputhp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    urputhp_cyl = np.zeros([nn_xx2, nn_r, nn_theta])
    for i in range(ti_n):
        ti = ti_arr[i]
        fn = 'DAT_{0:010d}.h5'.format(ti)
        print(fn)
        dict_vel = get_var_cart(fn, mode = 0)
        u_cart, v_cart, w_cart = dict_vel['u_cart'], dict_vel['v_cart'], dict_vel['w_cart']
        ux_cyl, ur_cyl, uth_cyl = get_vel_cyl(u_cart, v_cart, w_cart, zz*ll[2], yy, xx, cartpts_cyl)

        ux_cyl_m2 += ux_cyl
        ur_cyl_m2 += ur_cyl
        uth_cyl_m2 += uth_cyl
        uxpuxp_cyl += ux_cyl **2
        urpurp_cyl += ur_cyl **2
        uthputhp_cyl += uth_cyl **2
        uxpurp_cyl += ux_cyl * ur_cyl
        uxputhp_cyl += ux_cyl * uth_cyl
        urputhp_cyl += ur_cyl * uth_cyl    

    ux_cyl_m2 /= ti_n
    ur_cyl_m2 /= ti_n
    uth_cyl_m2 /= ti_n
    uxpuxp_cyl /= ti_n
    urpurp_cyl /= ti_n
    uthputhp_cyl /= ti_n
    uxpurp_cyl /= ti_n
    uxputhp_cyl /= ti_n
    urputhp_cyl /= ti_n

    uxpuxp_cyl -= ux_cyl_m2 ** 2
    urpurp_cyl -= ur_cyl_m2 ** 2
    uthputhp_cyl -= uth_cyl_m2 **2
    uxpurp_cyl -= ux_cyl_m2 * ur_cyl_m2
    uxputhp_cyl -= ux_cyl_m2 * uth_cyl_m2
    urputhp_cyl -= ur_cyl_m2 * uth_cyl_m2

    fn_out = 'umean_cyl.h5'
    f_out = h5py.File(fn_out, 'w')

    dset = f_out.create_dataset("ux", data = ux_cyl_m2)
    dset = f_out.create_dataset("ur", data = ur_cyl_m2)
    dset = f_out.create_dataset("uth", data = uth_cyl_m2)
    dset = f_out.create_dataset("uxpuxp", data = uxpuxp_cyl)
    dset = f_out.create_dataset("urpurp", data = urpurp_cyl)
    dset = f_out.create_dataset("uthputhp", data = uthputhp_cyl)
    dset = f_out.create_dataset("uxpurp", data = uxpurp_cyl)
    dset = f_out.create_dataset("uxputhp", data = uxputhp_cyl)
    dset = f_out.create_dataset("urputhp", data = urputhp_cyl)

    f_out.close()


# Mean momentum equation in cylindrical coordinate:
# \begin{equation}
# \frac{\partial u_x}{\partial t} + u_r \frac{\partial u_x}{\partial r} 
# + u_\theta \frac{\partial u_x}{\partial \theta} + u_x \frac{\partial u_x}{\partial x}
# = - \frac{\partial p}{\partial x} + \nu \frac{\partial^2 u_x}{\partial r^2}
# + \frac{\nu}{r}\frac{\partial u_x}{\partial r} + \frac{\nu}{r^2}\frac{\partial^2 u_x}{\partial \theta^2}
# + \nu \frac{\partial^2 u_x}{\partial x^2} + ReynoldsStress
# \end{equation}
# 
# \begin{equation}
# ReynoldsStress = - \frac{\partial \overline{u_x^\prime u_x^\prime}}{\partial x} 
# - \frac{1}{r}\frac{\partial (r \overline{u_x^\prime u_r^\prime})}{\partial r}
# - \frac{1}{r}\frac{\partial \overline{u_x^\prime u_\theta^\prime}}{\partial \theta}
# \end{equation}
# 
# Number the terms:
# 1. $u_r \frac{\partial u_x}{\partial r}$
# 2. $u_\theta \frac{\partial u_x}{\partial \theta}$
# 3. $u_x \frac{\partial u_x}{\partial x}$
# 4. $- \frac{\partial p}{\partial x}$
# 5. $\nu \frac{\partial^2 u_x}{\partial r^2}$
# 6. $\frac{\nu}{r}\frac{\partial u_x}{\partial r}$
# 7. $\frac{\nu}{r^2}\frac{\partial^2 u_x}{\partial \theta^2}$
# 8. $\nu \frac{\partial^2 u_x}{\partial x^2}$
# 9. $- \frac{\partial \overline{u_x^\prime u_x^\prime}}{\partial x} $
# 10. $- \frac{1}{r}\frac{\partial (r \overline{u_x^\prime u_r^\prime})}{\partial r}$
# 11. $- \frac{1}{r}\frac{\partial \overline{u_x^\prime u_\theta^\prime}}{\partial \theta}$

# In[13]:


i_post = False 
if i_post:
    fn_in = 'umean_cyl.h5'
    f_in = h5py.File(fn_in, 'r')    
    ux = np.array(f_in['ux'])
    ur = np.array(f_in['ur'])
    uth = np.array(f_in['uth'])
    uxpuxp = np.array(f_in['uxpuxp'])
    urpurp = np.array(f_in['urpurp'])
    uthputhp = np.array(f_in['uthputhp'])
    uxpurp = np.array(f_in['uxpurp'])
    uxputhp = np.array(f_in['uxputhp'])
    urputhp = np.array(f_in['urputhp'])
    f_in.close()
    
    # reference length, velocity, time
    #L0 = 1.0
    #U0 = 11.5258407161
    #T0 = L0 / U0
    #G0 = 9.81

    # Parameters for the compuational domain
    HBAR = 0.46
    #PEX = 1.45444104333
    #PEY = 8.72664625997

    # Parameters for wind profile
    kappa = 0.40
    Z0 = 3e-5
    RESBOT = 6987.45931461

    # Parameters for wave profile
    #HKA = 0.05
    #NWAVE = 12

    # Derived variables for wind
    RE = RESBOT * (1.0/kappa*np.log(HBAR / Z0))
    USBOT = 1.0 / (2.5 * np.log(HBAR/Z0))
    nu = 1.0 / RE
    
    fn_in = "mean_field_3d.h5"    
    f_in = h5py.File(fn_in, 'r') 
    p_cart = np.array(f_in['pp'])
    nut_cart = np.array(f_in['nut'])  
    f_in.close()
    interp_func = inp.RegularGridInterpolator((zz*ll[2], yy, xx), p_cart)
    p = interp_func(cartpts_cyl)
    interp_func = inp.RegularGridInterpolator((zz*ll[2], yy, xx), nut_cart)
    nut = interp_func(cartpts_cyl)
    nut += nu
    
    # term 1
    duxdr = grad_3d(ux, rr, axis = 1)
    term_1 = ur * duxdr

    # term 2
    duxdth = grad_3d(ux, thth, axis = 2)
    term_2 = uth * duxdth

    # term 3
    duxdx = grad_3d(ux, xx2, axis = 0)
    term_3 = ux * duxdx

    # term 4
    term_4 = - grad_3d(p, xx2, axis = 0)

    # term 5
    dduxdrdr = grad_3d(duxdr, rr, axis = 1)
    term_5 = nut * dduxdrdr

    # term 6
    term_6 = nut * duxdr
    for j in range(nn_r):
        term_6[:,j,:] /= rr[j]

    # term 7
    term_7 = nut * grad_3d(duxdth, thth, axis = 2)
    for j in range(nn_r):
        term_7[:,j,:] /= (rr[j]**2)

    # term 8
    term_8 = nut * grad_3d(duxdx, xx2, axis = 0)
    
    # term 9
    term_9 = - grad_3d(uxpuxp, xx2, axis = 0)
    
    # term 10
    temp = uxpurp.copy()
    for j in range(nn_r):
        temp[:,j,:] *= rr[j]
    term_10 = - grad_3d(temp, rr, axis = 1)
    for j in range(nn_r):
        term_10[:,j,:] /= rr[j]
        
    # term 11
    term_11 = - grad_3d(uxputhp, thth, axis = 2)
    for j in range(nn_r):
        term_11[:,j,:] /= rr[j]


