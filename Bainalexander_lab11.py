import sys 
import numpy as np
from matplotlib import pyplot as plt 


sys.path.append("C:\\Users\\Alex\\Documents\\GitHub\\BainAlexander_Lab10")


from BainAlexander_lab10 import make_tridiagonal
from BainAlexander_lab10 import make_initialcond
from BainAlexander_lab10 import spectral_radius




def advection1d(method, nspace, ntime, tau_rel,params):
    L, c = params
    x = np.linspace(-L / 2, L / 2, nspace)
    h = L/ nspace  
    tau_crit = h / c
    tau = tau_rel * tau_crit 
    t = np.arange(ntime) * tau 

    a = np.zeros((nspace,ntime))
    a[:, 0] = make_initialcond(sigma=0.2, K=35, x=x)
    for n in range(0,ntime-1):
        if method == "FTCS":
            A = make_tridiagonal(nspace, 1 - tau/h, 1, -1 + tau/h)
        elif method == "Lax":
            A = make_tridiagonal(nspace, 0.5 - tau/(2*h), 0, 0.5 + tau/(2*h))
        else: 
            return "INVALID METHOD"
    
    if spectral_radius(A) > 1:
        print("Warning: Unstable integration requested!")

    for n in range(1,ntime):
        a[:, n+1] = np.dot(A, a[:, n])
    return a,x,t
 
nspace = 300
ntime = 500
tau_rel = 1
params = [5, 1]  # L=5, c=1

a, x, t = advection1d("Lax", nspace, ntime, tau_rel, params)


plt.figure(figsize=(10, 6))
for n in range(0, ntime, ntime // 10):  
    plt.plot(x, a[:, n], label=f'Time = {t[n]:.2f}')
plt.xlabel('x')
plt.ylabel('Wave Amplitude a(x, t)')
plt.title('Wave Propagation (Lax Method)')
plt.legend()
plt.grid()
plt.show()