import sys 
import numpy as np
from matplotlib import pyplot as plt 

#this makes python check in the right spot for my lab10.py code 
#because it is in a different folder 
sys.path.append("C:\\Users\\Alex\\Documents\\GitHub\\BainAlexander_Lab10")

#imported functions from previous lab
from BainAlexander_lab10 import make_tridiagonal
from BainAlexander_lab10 import make_initialcond
from BainAlexander_lab10 import spectral_radius



#this is the advection1d function we created for part 1 
def advection1d(method, nspace, ntime, tau_rel,params):
    #these are our parameters
    L, c = params
    x = np.linspace(-L / 2, L / 2, nspace)
    h = L/ nspace  
    tau_crit = h / c
    tau = tau_rel * tau_crit 
    t = np.arange(ntime) * tau 
    #creates array and adds initial conditions to it based on function
    #created in lab 10 
    a = np.zeros((nspace,ntime))
    a[:, 0] = make_initialcond(sigma=0.2, K=35, x=x)
    #for loop that runs for each timestep 
    for n in range(0,ntime-1):
        #this is the calculation for FTCS method
        if method == "FTCS":
            A = make_tridiagonal(nspace, 1 - tau/h, 1, -1 + tau/h)
            #calculation for Lax method
        elif method == "Lax":
            A = make_tridiagonal(nspace, 0.5 - tau/(2*h), 0, 0.5 + tau/(2*h))
            #returns an error if an invalid method is added 
        else: 
            return "INVALID METHOD"
    #this creates a warning for instability based off the spectral radius function
    #we created in lab10
    if spectral_radius(A) > 1:
        print("Warning: Unstable integration requested!")
    #iterates over the time steps except the last one 
    for n in range(0,ntime-1):
        #this updates the scalar field 
        a[:, n+1] = np.dot(A, a[:, n])
    return a,x,t
 
#these are some values we give to the function
nspace = 300
ntime = 500
tau_rel = 1
params = [5, 1]  # L=5, c=1

a, x, t = advection1d("Lax", nspace, ntime, tau_rel, params)

#here we plot the info 
plt.figure(figsize=(10, 6))
#for loop that plots the wave at 10 evenly spaced time steps
for n in range(0, ntime, ntime // 10):  
    plt.plot(x, a[:, n], label=f'Time = {t[n]:.2f}')
#these are axis labels 
plt.xlabel('x')
plt.ylabel('Wave Amplitude a(x, t)')
plt.title('Wave Propagation (Lax Method)')
plt.legend()
plt.grid()
plt.show()
