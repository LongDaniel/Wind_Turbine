import numpy as np
import os

fn = '9030to33000'
fn1 = '../'+fn
cmd = 'mkdir -p '+fn1
print(cmd)
os.system(cmd)

cmd = 'mv -t '+fn1+' mean_* POST_UM_* log.post.* spectrum_mean.dat *_tm.h5'
print(cmd)
os.system(cmd)

for i in range(1,16):
    fn2 = 'POST_U_1D1_{0:04d}'.format(i)
    cmd = 'mkdir -p '+fn1+'/'+fn2
    print(cmd)
    os.system(cmd)
    
    cmd = 'cp -t '+fn1+'/'+fn2+' '+fn2+'/post_u_1d_s1.log'+' '+fn2+'/result_mean.npz'
    print(cmd)
    os.system(cmd)

