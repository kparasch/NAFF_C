from __future__ import print_function
import sys
sys.path.append('..')
sys.path.append('.')
import nafflib as NAFF

import numpy as np

np.random.seed(123456)
N=100
noise_rms=0.1
i = np.linspace(1,N,N)

q_true=6.24783351

x=np.empty_like(i,dtype=np.complex128)
print('                  True     frequency is Q_true = {0:.10f}'.format(q_true))

x  = 1*np.cos(2*np.pi*q_true*i)*np.exp(-i*2.3) 
tune, tau = NAFF.get_tune_tau(x);
print('(real           ) Estimated frequency is Q_hat = {0:.10f}'.format(tune))
print('(real           ) Estimated damping is tau_hat = {0:.10f}'.format(tau))

