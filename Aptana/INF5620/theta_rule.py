from dolfin import *
import numpy as np
from math import log
def diffusion(dt,size):
    
    # Create mesh and define function space
    #mesh = UnitSquareMesh(size, size)
    mesh = UnitIntervalMesh(size)
    V = FunctionSpace(mesh, "Lagrange", 1)
    R = FunctionSpace(mesh, "R", 0)
    W = V * R
    
    # Define variational problem
    (u, c) = TrialFunction(W)
    (v, d) = TestFunctions(W)
    #f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
    #g = Expression("-sin(5*x[0])")
    f = Constant(0.0)
    g = Constant(0.0)
    #dt = 0.001
    T = 0.05
    rho = 1.
    I = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0)

    u_1 = interpolate(I,V) 
    alpha = Constant(1.0)
    #a = ((dt/rho)*inner(alpha*grad(u), grad(v)) + (dt/rho)*c*v + u*d)*dx +u*v*dx
    #L = f*v*dx + g*v*ds + u_1*v*dx
    F = (rho/dt)*(u-u_1)
    
    # Compute solution
    w = Function(W)
    t = dt
    while t<=T:
        I.t = t
        solve(a == L, w)
        (u, c) = w.split(True)
        u_1.assign(u)
        t+=dt
    
    #plot(u)
    
    u_e = interpolate(I, V)
    #plot(u_e)
    #plot(u,interactive=True)
    #plot(u_e,interactive=True)
    e = np.abs(u_e.vector().array()-u.vector().array()).max()
        
    E = np.sqrt(np.sum(e**2)/u.vector().array().size)
    return e

# h = [0.1,0.01,0.001]
# size = [5,50,500]
# #dt = 0.0001
# #size = 5000
# E = []
# Eh = []
# for dt,s in zip(h,size):
#     e = diffusion(dt,s)
#     E.append(e)
# for i,j in zip(E,h):
#     Eh.append(i/j)
# print Eh
# print E
h = [0.00045,0.00040,0.00035,0.00030]
h_hf = [k/2 for k in h]
 
# error due to h
size = 64
e = []
for dt in h:
    E = diffusion(dt,size)
    e.append(E)
# error du to h/2
e_hf = []
for dt in h_hf:
    E = diffusion(dt,size)
    e_hf.append(dt)
     
P = []
Eh = []
for i,j in zip(h,h_hf):
    p = log(i/j)/log(2)
    P.append(p)
 
# change dt and dx
ee = []
size = [3100,2200,1300,1040]
#size = [4000,3000,2000,1000]
for dt,s in zip(h,size):
    E = diffusion(dt,s)
    ee.append(E)
 
for i,j in zip(ee,h):
    eh = i/j
    Eh.append(eh)
     
print Eh 
print P
 
print " eror", ee
#     


    
    





