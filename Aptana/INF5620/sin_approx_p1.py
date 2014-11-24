

from math import pi, sin, cos
import numpy as np

import matplotlib.pyplot as pl
from scipy.integrate import quad


   
def sin_approx_P1():
    """
    Approximate sin(x) by P1 finite elements.
    We use the reference basis function:
    phi0(X) = (1/2)*(1-X)
    phi1(X) = (1/2)*(1+X)
    """
    
    # global matrixes A and b
    N = 3
    h = pi/2
    xm0 =  (0+pi/2)/2
    xm1 =  (pi/2. + pi)/2. 
    A = np.zeros((N,N))
    b = np.zeros(N)
    
    # local matrixes Ae and be
    Ae = np.zeros((2,2))
    be = np.zeros(2)
    bei = np.zeros(2)
    Ae[0,0], Ae[0,1], Ae[1,1], Ae[1,0] = h/3., h/6, h/3, h/6
    
    
    def f(X):
        h = pi/2.
        xm0 = (0+pi/2)/2
        return sin(xm0+0.5*h*X)*(1-X)*(h/4)
    def g(X):
        h = pi/2.
        xm1 =  (pi/2. + pi)/2. 
        return sin(xm1+0.5*h*X)*(1+X)*(h/4)
        
    bei[0], err = quad(f,-1,1)
    bei[1], err = quad(g,-1,1)
    
    be[0] = cos(xm0-0.5*h)-(1/h)*(sin(xm0+0.5*h)-sin(xm0-0.5*h) )
    be[1] = -cos(xm1+0.5*h)+(1/h)*(sin(xm1+0.5*h)-sin(xm1-0.5*h) )
    
    #### assemble element matrices###
    elements = [ [0,1],[1,2] ]
    
    # number of elements
    E = len(elements)
    coordinates = [i*(pi/2) for i in range(E+1)]
    
    Id = [i for i in range(E)]
    
    for e in range(E):
        for r in Id:
            for s in Id:
                i = elements[e][r]
                j = elements[e][s]
                
                A[i,j] += Ae[r,s]
                b[i]   += be[r]
                
    c = np.linalg.solve(A,b)
    
    print c
    
    
    
    

sin_approx_P1()
    
    
    