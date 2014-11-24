"""
solve the time dependent poisson equation with 
Newmann boundary condition:
rho*u_t = div(alpha(u)grad(u)) +f
du/dn = 0.
Because of the Newmann boundary condition the solution u is known to a constant c.
To find this constant extra condition must be apply:
integral(u)dx = integral(du/dn)*ds. 

Now we  must solve the problem:
rho*u_t = div(alpha(u)grad(u))+c +f

subject to the constraint:
integral(u)dx = integral(du/dn)*ds. 

"""

from dolfin import *
import numpy as np
from math import log
def solver(dimension,d,I,f,rho,size,T,dt,alpha=None):
    """
    usage of this function:
    
    solver(dimension,d,alpha,I,f,rho,size) return the 
    solution of the time dependent poisson equation.
    dimension is the dimension of the problem
    d is the degree of Lagrange element: d= 1 or d = 2
    alpha is the nonlinear term in the equation
    I is the initial condition
    f is the source term 
    rho is the density
    size is the mesh size
    """
    
    #specify dimension of the problem
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    
    elif  dimension == 2:
        mesh = UnitSquareMesh(size,size)
        
    elif  dimension == 3:
        mesh = UnitCubeMesh(size,size,size)
        
    # mixt finite element function space formulation
    V = FunctionSpace(mesh,'Lagrange',d) 
    R = FunctionSpace(mesh, "R", 0) # space for real number
    W = V * R                       # cross product of V and R
    
    # Define variational problem
    (u, c) = TrialFunction(W)
    (v, d) = TestFunctions(W)
    g = Constant(0.0)
    #dt = 0.001
    #T = 0.05
    #rho = 1.
    u_1 = interpolate(I,V) 
    if alpha == None:
        alpha = interpolate(I,V)
    a = ((dt/rho)*inner(alpha*grad(u), grad(v)) + (dt/rho)*c*v + u*d)*dx +u*v*dx
    L = f*v*dx + g*v*ds + u_1*v*dx
    
    # Compute solution
    w = Function(W)
    t = dt
    while t<=T:
        I.t = t
        solve(a == L, w)
        (u, c) = w.split(True)
        u_1.assign(u)
        t+=dt
    
    u_e = interpolate(I, V) # analytical solution
    e = np.abs(u_e.vector().array()-u.vector().array()).max() # max difference 
    L2_norm = np.sqrt(np.sum(e**2)/u.vector().array().size)         # L2 norm
    
    return u,u_e,e,L2_norm
    


 
def constant_solution_verification():
    C = 0
    I = Constant(C)
    f = Constant(0.0)
    rho = 1.0
    d = 2
    dimension = 1
    T = 0.05
    dt = 0.0003
    size = 1040
    #alpha = Constant(1.0)
    u_num,u_anal,error,L2_norm = solver(dimension,d,I,f,rho,size,T,dt)
    
    tol = 1e-15
    assert error <= tol
    print error
    
def analytical_solution_verification():
    I = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0)
    f = Constant(0.0)
    rho = 1.
    d = 1
    T = 0.05
    dt = 0.0003
    size = 1040
    alpha = Constant(1.0)
    dimension = 1
    u_num,u_anal,error,L2_norm = solver(dimension,d,I,f,rho,size,T,dt,alpha)
    tol = 1e-4
    assert L2_norm <=tol
    

def convergence_analytical_solution():
    I = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0)
    f = Constant(0.0)
    rho = 1.
    d = 1
    T = 0.05
    #dt = 0.0003
    #size = 1040
    alpha = Constant(1.0)
    dimension = 1
    
    h = [0.00045,0.00040,0.00035,0.00030] #dt list
    h_hf = [k/2 for k in h] # dt/2 list
     
    # error due to h
    size = 64
    e = []
    for dt in h:
        u_num,u_anal,E,L2_norm = solver(dimension,d,I,f,rho,size,T,dt,alpha)

        e.append(E)
    # error du to h/2
    e_hf = []
    for dt in h_hf:
        u_num,u_anal,E,L2_norm = solver(dimension,d,I,f,rho,size,T,dt,alpha)
        e_hf.append(dt)
         
    P = []
    Eh = []
    for i,j in zip(h,h_hf):
        p = log(i/j)/log(2)
        P.append(p)
     
    # change dt and dx
    ee = []
    siz = [3100,2200,1300,1040]
    #size = [4000,3000,2000,1000]
    for dt,s in zip(h,siz):
        u_num,u_anal,E,L2_norm = solver(dimension,d,I,f,rho,s,T,dt,alpha)
        ee.append(E)
     
    for i,j in zip(ee,h):
        eh = i/j
        Eh.append(eh)
         
#     print Eh 
#     print P
#      
#     print " eror", ee
    #     
        
  
#convergence_analytical_solution()        
constant_solution_verification()
#analytical_solution_verification()
        
        
        
        
        
        