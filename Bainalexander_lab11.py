import sys 
import numpy as np
from matplotlib import pyplot as plt 


sys.path.append("C:\\Users\\Alex\\Documents\\GitHub\\BainAlexander_Lab10")

from BainAlexander_lab10 import make_tridiagonal, make_initialcond

a = make_tridiagonal(5,3,1,4)
print(a)

def advection1d(method, nspace, ntime, tau_rel,params):
    L, c = params
    x = np.linspace(-L / 2, L / 2, nspace)
    h = x[1] - x[0]  
    tau_crit = h / abs(c)  
    tau = tau_rel * tau_crit 
    t = np.arange(ntime) * tau 