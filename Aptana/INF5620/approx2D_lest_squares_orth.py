"""
approximation of the function:
f(x,y)= x*(1-x)*y*(1-y)*exp(-x-y) by using an orthonormal basis function:
phi_i(x,y) = sin(pi*p*x)*sin(q*pi*y). phi_i(x,y) is an orthonormal basis function if

(phi_i,phi_j) = 0 for i != j.
"""
import pylab as pl
from pylab import *
import matplotlib.pyplot as plt
from scipy import integrate
from mpl_toolkits.mplot3d.axes3d import Axes3D
import numpy as np
import sympy as sp
def least_square_orth(Nx,Ny,f,omega,phi,string):
    """
    return the lest square approximation of a function f in 2D.
    with an orthonormal basis function. The function also return 
    the maximum difference between f and it least square appromimation.
    
    usage: f,error = least_square_orth(Nx,Ny,f,omega,phi)
    
    out put:
    ######################################################################
    Nx =  3 ; Ny =  3
    finite element approximation at x, y  = 0.5 , 0.1:    0.0116711446489
    exact solution at x,y = 0.5, 0.1                :    0.0123482618121
    ######################################################################

    """
    x, y = sp.symbols("x y")
    # define index arrays for p and q
    I_x, I_y = [i for i in range(Nx+1)], [i for i in range(Ny+1)]
    
    # define index array for i = q*(Nx+1)+p:
    index_i = [q*(Nx+1)+p for q in I_y for p in I_x]

    # define matrix A and b
    N = len(index_i) # size of A and b
    A, b = np.zeros((N,N)), np.zeros(N)
    
    if string == "orthogonal":
        for i in range(N):
            integrand_phi_f  = phi[i]*f
            integrand2   = sp.lambdify([x,y],integrand_phi_f)
            b[i] = sp.mpmath.quad(integrand2,[omega[0][0],omega[0][1]], [omega[1][0],omega[1][1]] )
            integrand_phi_phi = phi[i]*phi[i]
            integrand1 = sp.lambdify([x,y],integrand_phi_phi) 
            A[i,i] = sp.mpmath.quad( integrand1,[omega[0][0],omega[0][1]], [omega[1][0],omega[1][1]]   )
    else:
        
        for i in range(N):
            integrand_phi_f  = phi[i]*f
            integrand2   = sp.lambdify([x,y],integrand_phi_f)
            b[i] = sp.mpmath.quad(integrand2,[omega[0][0],omega[0][1]], [omega[1][0],omega[1][1]] )
            for j in range(N):
                integrand_phi_phi = phi[i]*phi[j]
                integrand1 = sp.lambdify([x,y],integrand_phi_phi)
                
                A[i,j] = sp.mpmath.quad( integrand1,[omega[0][0],omega[0][1]], [omega[1][0],omega[1][1]]   )
    
        
            

    #c = np.linalg.solve(A,b)
    c = [b[i]/A[i,i] for i in range(len(b))]
    g = sum([c[i]*phi[i] for i in index_i])
    
    # approximation of function f
    fapprox = sp.lambdify([x,y],g,"numpy")
    fe = sp.lambdify([x,y],f,"numpy")
    
    # error between f and it approximation
    X = np.linspace(0,1,100)
    Y = np.linspace(0,1,100)
    print"######################################################################"
    print"Nx = ",Nx,";", "Ny = ",Ny
    print"finite element approximation at x, y  = 0.5 , 0.1:   ", fapprox(0.5,0.1)
    print "exact solution at x,y = 0.5, 0.1                :   ", fe(0.5,0.1) 
    print"######################################################################"
    X1, Y1 = pl.meshgrid(X,Y)
    U = fapprox(X1,Y1)
    U1 = fe(X1,Y1)
    fig = pl.figure(figsize=(14,6))
    ax = fig.add_subplot(1,2,1, projection='3d')
    ax.plot_surface(X1, Y1, U, rstride=1, cstride=1, alpha=0.5, cmap= cm.coolwarm)
    pl.title("finite element approximation")
    ax1 = fig.add_subplot(1,2,2, projection='3d')
    ax1.plot_surface(X1, Y1, U1, rstride=1, cstride=1, alpha=0.5,cmap=cm.coolwarm, linewidth=0, antialiased=False)
    pl.title("exact solution")
    pl.show()
    

omega = [[0,1],[0,1]]
Nx, Ny = 3,3
x,y = sp.symbols("x y")

I = [q*(Nx+1)+p for q in range(Ny+1) for p in range(Nx+1)]

phi = [sp.sin(sp.pi*(p+1)*x)*sp.sin(sp.pi*(q+1)*y)  for q in range(Ny+1) for p in range(Nx+1)]

#phi = [(x**p)*(y**q) for p in range(Ny+1) for q in range(Nx+1)]
f = x*(1-x)*y*(1-y)*sp.exp(-x-y)
#f = (1+x**2)*(1+2*y**2)

f = least_square_orth(Nx,Ny,f,omega,phi,"orthogonal")


    
    
    
    
    
    
    