
import numpy as np
from math import pi, tanh, sin
import matplotlib.pyplot as plt

from scipy.integrate import quad
import sympy as sp
import sympy.mpmath  as sm



def approximation(k):
    """
    Approximate sin(x) by 
    1) using Galerkin or projection method f. The projection method 
    gives the same result as the least square method.
    """
    N = 8
    q = 100 # number mesh point for plot
    A = np.zeros((N+1,N+1)) 
    b = np.zeros(N+1)
    
    
    t = sp.symbols("t")      

    for i in range(N+1):
        phiiAndf = lambda x: (x**(2*i+1))*sin(x)
        b[i],err =quad( phiiAndf,0,k*pi*0.5  ) 
        
        for j in range(N+1):
            phiiAndphij = lambda x: x**(2*i+1)*x**(2*j+1)
            #b[i] = sp.integrate(phi[i]*f, (x,0,2*pi) )
           
            A[i,j],err = quad(phiiAndphij, 0,k*pi*0.5  )
            
    
    #solve Ac = b for c
    c = np.linalg.solve(A, b)
    
    # approximation of sin(t) by projection method
    phi = [t**(2*i+1) for i in range(N+1)]
    u = sum([ c[i]*phi[i] for i in range(N+1)])
    h = sp.lambdify(t,u,"numpy")
    z = np.linspace(0,k*pi*0.5,q)
    
    # Taylor series
    def g(r):
        return sin(r)
    
    
    z = np.linspace(0,k*pi*0.5,q)
    z0 = np.linspace(0,k*pi*0.5,q)
    deg = 9
    
    Taylor = []
    for x,x0 in zip(z,z0):
        p = sm.taylor(g,x0,deg)
        r = sm.polyval(p[::-1], x-x0)
        Taylor.append(r)   
    
    plt.plot(z,h(z),"ro",z,Taylor,"b--", z,np.sin(z))
    #plt.plot(z,Taylor)
    plt.show()

    
if __name__=="__main__":
    approximation(4)
    
    
    
    
    
    
    
    