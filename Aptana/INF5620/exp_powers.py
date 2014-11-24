import numpy as np
import pylab as pl
from sympy import *
import sys
import math
import sympy.mpmath  as sm


def approximate(N):
    """
    Approximate exp(-x) by 
    1) using Galerkin or projection method for N = 2,4,6. The projection method 
    gives the same result as the least square method.
    
    2) using Taylor polynomial of degree N
    for N = 2,4,6
    plot the different approximation on the same figure
    
    """
    #c = np.zeros(N+1)  # linear combination coefficient
    A = np.zeros((N+1,N+1)) 
    b = np.zeros(N+1)
    
    
    x = Symbol("x")
    
    # basis functions
    phi = [x**i for i in range(N+1)] 
    
    # matrix A and b
    for i in range(N+1):
        for j in range(N+1):
            b[i] = integrate(phi[i]*(exp(-x)), (x,0,4))
            A[i,j] = integrate(phi[i]*phi[j], (x,0,4))
            
    
    #solve Ac = b for c
    c = np.linalg.solve(A, b)
     
    # approximation of exp(-x) by projection method
    u = sum([ c[i]*phi[i] for i in range(N+1)])
    f = lambdify(x,u)
    
    #taylor polynomial approximation of exp(-x)
    def g(r):
        return exp(-r)
    
    q = 100
    z = np.linspace(0,4,q)
    x0  = 2
    deg = N
    p = sm.taylor(g,x0,deg)
    
    Taylor = []
    for x in z:
        r = sm.polyval(p[::-1], x-x0)
        Taylor.append(r)   
     
    # ploting 
    pl.plot(z,f(z),"ro",z,np.exp(-z), z,Taylor,"g--")
    pl.legend([" Galerkin N = %g" %N, "exp(-x)","Taylor "])
    pl.show()
#         
if __name__=="__main__":
    approximate(6)
    
    
    
    