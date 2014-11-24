from dolfin import *
import numpy as np

def diffusion(I,f,dt,T,t,rho,size,degree,d,dimension,alpha = None):
    
    # Create mesh and define function space
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    elif dimension == 2:
        mesh == UnitSquareMesh(size,size)
    elif dimension == 3:
        mesh = UnitCube(size,size,size)
    
    # function space
    V = FunctionSpace(mesh, "Lagrange", d)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_1 = interpolate(I,V) 
    if alpha == None:
        
        def alpha(u):
            return 1
    
    #the variational formulation for Picard method
    a = u*v*dx +dt/rho*inner(alpha(u_1)*grad(u),grad(v))*dx
    L = (u_1 + f*dt/rho)*v*dx
    t = dt
    while t<=T:    
        #update f and solve the variational problem
        f.t = t
        u = Function(V)
        solve(a == L, u)
        #update values for next iteration
        t += dt
        u_1.assign(u)

    return  u
    