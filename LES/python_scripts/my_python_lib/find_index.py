import numpy as np

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

if __name__ == "__main__":
	z_test = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
	limits_test = [0.23, 0.61]
	imin_test, imax_test = find_index (z_test, limits_test)
	print('z = '+str(z_test))
	print('limits = '+str(limits_test))
	print('imin, imax = '+str([imin_test, imax_test]))
	print('values inside the limits are '+str(z_test[imin_test:imax_test]))