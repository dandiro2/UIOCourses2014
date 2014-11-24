from dolfin import *
import numpy as np
import pylab as pl

def diffusion(dimension):
    
    # Create mesh and define function space
    size = 120
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    elif dimension == 2:
        mesh = UnitSquareMesh(size,size)
    elif dimension == 3:
        mesh = UnitCube(size,size,size)
    

    rho = 1.0
    d = 1
    
    def alpha(u):
        return 1 + u**2
        #alpha = lambda u: 1 + u**2
    f = Expression('-rho*pow(x[0],3)/3 + rho*pow(x[0],2)/2 \
    + 8*pow(t,3)*pow(x[0],7)/9 - 28*pow(t,3)*pow(x[0],6)/9 \
    + 7*pow(t,3)*pow(x[0],5)/2 - 5*pow(t,3)*pow(x[0],4)/4 \
    + 2*t*x[0] - t', rho=rho, t=0)
    I = Expression('exp(-(1./2*sigma*sigma)*(x[0]*x[0]+x[1]*x[1]) )',sigma = 1.0)
    V = FunctionSpace(mesh, "Lagrange", 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_1 = interpolate(I,V) 
    
    T = 3
    dt = 0.01
    #the variational formulation for Picard method
    a = u*v*dx +dt/rho*inner(alpha(u_1)*grad(u),grad(v))*dx
    L = (u_1 + f*dt/rho)*v*dx
    t = dt
    while t<=T:    
        #update f and solve the variational problem
        f.t=t
        u = Function(V)
        solve(a == L, u)
        #update values for next iteration
        t += dt
        u_1.assign(u)

    
    plot(u)
    interactive()

dimension = 2
diffusion(dimension)
