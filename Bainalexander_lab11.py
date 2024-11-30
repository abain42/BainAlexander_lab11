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
    h = L/ nspace  
    tau_crit = h / c
    tau = tau_rel * tau_crit 
    t = np.arange(ntime) * tau 

    for n in range(0,ntime-1):
        if method == "FTCS":
            A = make_tridiagonal(nspace, 1 - tau/h, 1, -1 + tau/h)
        elif method == "Lax":
            A = make_tridiagonal(nspace, 0.5 - tau/(2*h), 0, 0.5 + tau/(2*h))
        else: 
            return "INVALID METHOD"