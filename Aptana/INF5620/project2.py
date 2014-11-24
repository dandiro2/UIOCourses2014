"""
solving the nonlinear diffusion equation
with homogeneaous Newmann boundary condition
"""


from dolfin import *
import numpy as np
from math import pi,cos,exp
def diffusion(dimension,d,alpha,I,f,rho,size):
    """
    solve the time dependent diffusion equation
    with homogeneous Newmann boundary condition.
    in the program u_1 correspond to u^k and u correspond
    to u^k+1.
    usage of the function duffusion:
    
    diffusion('dimension','d',I,alpha,rho):
    dimension = 1D, 2D, 3D for a one dimension,
    two dimension, and three dimension. problem respectivelly 
    
    d = 1 for P1 Lagrange element and d = 2 for P2 Lagrange 
    element.
    """
    
    # specify the dimension of the problem
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    
    elif  dimension == 2:
        mesh = UnitSquareMesh(size,size)
        
    elif  dimension == 3:
        mesh = UnitCubeMesh(size,size,size)
        
    
    # finite element function space
    V = FunctionSpace(mesh,'Lagrange',d)
    
    u_1 = interpolate(I,V)
    u = TrialFunction(V)
    v = TestFunction(V)
    
    dt = 0.02
    a = u*v*dx + (dt/rho)*inner(alpha*grad(u),grad(v))*dx 

    L = u_1*v*dx + (dt/rho)*f*v*dx
    A = assemble(a)
    b = assemble(L)
    
    #Compute solution
    u= Function(V)
    T = 0.05
    t = dt
    while t <= T:
        I.t = t
        solve(A,u.vector(),b)
        t += dt
        u_1.assign(u)
    return u



def verification_constant(dimension,d,size):
    C = 0.5
    f = Constant(0.0)
    rho = 1.
    
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    
    elif  dimension == 2:
        mesh = UnitSquareMesh(size,size)
        
    elif  dimension == 3:
        mesh = UnitCubeMesh(size,size,size)
    V = FunctionSpace(mesh,'Lagrange',d)
    u_ex = Constant(C)
    u_exact = interpolate(u_ex,V)
    I = Constant(C)
    #set alpha value at = 0
    alpha = interpolate(I,V)
    U = diffusion(dimension,d,alpha,I,f,rho,size)
    er = np.abs(U.vector().array()-u_exact.vector().array()).max()
    E = np.sqrt(np.sum(er**2)/U.vector().array().size)
    #print U.vector().array()
    print er
    print "L2 norm",E
    #plot(U,interactive=True)

def verification_analytical(dimension,d,size):
    
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    
    elif  dimension == 2:
        mesh = UnitSquareMesh(size,size)
        
    elif  dimension == 3:
        mesh = UnitCubeMesh(size,size,size)
    V = FunctionSpace(mesh,'Lagrange',d)
    alpha = Constant(1.0)
    f = Constant(0.0)
    rho = 1
    I = Expression('cos(pi*x[0])')
    u_ex = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0.05)
    u_exact = interpolate(u_ex,V)
    
    U = diffusion(dimension,d,alpha,I,f,rho,size)
    
    e = U.vector().array()-u_exact.vector().array()
    
    E = np.sqrt(np.sum(e**2)/U.vector().array().size)
    print E
    print abs(e).max()
    plot(u_exact)
    plot(U,interactive=True)


  
#verification_constant(dimension=1,d=1,size=100)
verification_analytical(dimension=1,d=2,size=100)

