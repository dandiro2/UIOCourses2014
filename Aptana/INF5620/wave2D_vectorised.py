"""
2a82a0785570545692c68a4d6c
solve the 2D wave equation with variable wave velocity q(x,y).
d2u/dt2 + b(du/dt) = d/dx(q(x,y)du/dx) + d/dy(q(x,y)du/dy) +f(x,y,t)
BC: du/dn = 0 (Newman BC)
IC: du(x,y,0)/dt = V(x,y)
    u(x,y,0)     = I(x,y)
"""
from pylab import *
import numpy as np
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D
import sympy as sp
from math import pi
def vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b):
    """"
    implementation of the vectorized solver.
    I(x,y) is the initial solution, V(x,y) is the initial velocity (at t=0), q(x,y) is the 
    wave velocity at subsequent time. Lx and Ly are the length of the domain along the x and y
    axis. Nx and Ny are number of points along the x and y axis. dt is the time step, T is the 
    total simulation time. ste1BC and next_stepBC are the boundary condition at the first time
    step and the subsequent time step respectively.
    """
    
    u   = np.zeros((Nx+1,Ny+1)) # solution array
    u_1 = np.zeros((Nx+1,Ny+1)) # solution u[n]
    u_2 = np.zeros((Nx+1,Ny+1)) # solution u[n-2]
    f_ = np.zeros((Nx+1,Ny+1))
    V_ = np.zeros((Nx+1,Ny+1))
    x  = np.linspace(0, Lx, Nx+1) # mesh points in x dir
    y  = np.linspace(0, Ly, Ny+1) # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    xv, yv = x[:,np.newaxis],y[np.newaxis,:] # used to vectorise f, V and I


    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1) # mesh points in time
    
    # simplification parameters k1,k2,k3,r1,r2,r3:
    k1 = ( (dt*b/2.)-1 )/( (dt*b/2.)+1 )
    k2 = 2./(  (dt*b/2.)+1   )
    k3 = dt**2/((dt*b/2.)+1)
    
    r1 = np.zeros((Nx+1,Ny+1))
    #V  = np.zeros((Nx+1,Ny+1))
    q  = np.zeros((Nx+1,Ny+1))
    
    q[0:Nx+1,0:Ny+1] = velocity(xv,yv)
    r2 = k2/(1-k1)
    r3 = k3/(1-k1)
    
    import time;  t0 = time.clock()  # for measuring CPU time
    
    # initial condition
    u_1[:,:] = I(xv,yv)
    
    
    #special formula for the first time step
    n = 0
    
    rx,ry = (1./dx**2),(1./dy**2)
    Dxx = (0.5*(q[1:Nx,1:Ny]+q[2:Nx+1,1:Ny])*(u_1[2:Nx+1,1:Ny]-u_1[1:Nx,1:Ny])-0.5*(q[1:Nx,1:Ny]+q[0:Nx-1,1:Ny])*(u_1[1:Nx,1:Ny]-u_1[0:Nx-1,1:Ny]))*rx
    Dyy = (0.5*(q[1:Nx,1:Ny]+q[1:Nx,2:Ny+1])*(u_1[1:Nx,2:Ny+1]-u_1[1:Nx,1:Ny])-0.5*(q[1:Nx,1:Ny]+q[1:Nx,0:Ny-1])*(u_1[1:Nx,1:Ny]-u_1[1:Nx,0:Ny-1]))*ry
    delta = Dxx+Dyy
    f_[:,:] = f(xv,yv,t[n])
    V_[:,:] = V(xv,yv)
    u[1:Nx,1:Ny]=( (2*k1*dt)/(k1-1)  )*V_[1:Nx,1:Ny]+r2*u_1[1:Nx,1:Ny]+r3*(delta + f_[1:Nx,1:Ny] )
    
            
      
    # boundary condition
       
    #i = 0 and 0<= j < Ny
    Dxx = (0.5*(q[0,0:Ny]+q[1,0:Ny])*(u_1[1,0:Ny]-u_1[0,0:Ny])-0.5*(q[0,0:Ny]+q[1,0:Ny])*(u_1[0,0:Ny]-u_1[1,0:Ny]))*(1./dx**2)
    Dyy = (0.5*(q[0,0:Ny]+q[0,1:Ny+1])*(u_1[0,1:Ny+1]-u_1[0,0:Ny])-0.5*(q[0,0:Ny]+q[0,1:Ny+1])*(u_1[0,0:Ny]-u_1[0,1:Ny+1]))*(1./dy**2)
    delta = Dxx+Dyy
    u[0,0:Ny] = ( (2*k1*dt)/(k1-1)  )*V_[0,0:Ny]+r2*u_1[0,0:Ny]+r3*(delta + f_[0,0:Ny] )
    # i = 0 and  j = Ny
    Dxx = (0.5*(q[0,Ny]+q[1,Ny])*(u_1[1,Ny]-u_1[0,Ny])-0.5*(q[0,Ny]+q[1,Ny])*(u_1[0,Ny]-u_1[1,Ny]))*(1./dx**2)
    Dyy = (0.5*(q[0,Ny]+q[0,Ny-1])*(u_1[0,Ny-1]-u_1[0,Ny])-0.5*(q[0,Ny]+q[0,Ny-1])*(u_1[0,Ny]-u_1[0,Ny-1]))*(1./dy**2)
    delta = Dxx+Dyy
    u[0,Ny] = ( (2*k1*dt)/(k1-1)  )*V_[0,Ny]+r2*u_1[0,Ny]+r3*(delta +f_[0,Ny] )
    
    
    #j = 0 and 0<= i < Ny
    Dxx = (0.5*(q[0:Nx,0]+q[0:Nx,1])*(u_1[0:Nx,1]-u_1[0:Nx,0])-0.5*(q[0:Nx,0]+q[0:Nx,1])*(u_1[0:Nx,0]-u_1[0:Nx,1]))*(1./dx**2)
    Dyy = (0.5*(q[0:Nx,0]+q[1:Nx+1,0])*(u_1[1:Nx+1,0]-u_1[0:Nx,0])-0.5*(q[0:Nx,0]+q[1:Nx+1,0])*(u_1[0:Nx,0]-u_1[1:Nx+1,0]))*(1./dy**2)
    delta = Dxx+Dyy
    u[0:Nx,0] = ( (2*k1*dt)/(k1-1)  )*V_[0:Nx,0]+r2*u_1[0:Nx,0]+r3*(delta + f_[0:Nx,0])
    # j = 0 and  i = Nx
    Dxx = (0.5*(q[Nx,0]+q[Nx,1])*(u_1[Nx,1]-u_1[Nx,0])-0.5*(q[Nx,0]+q[Nx,1])*(u_1[Nx,0]-u_1[Nx,1]))*(1./dx**2)
    Dyy = (0.5*(q[Nx,0]+q[Nx-1,0])*(u_1[Nx-1,0]-u_1[Nx,0])-0.5*(q[Nx,0]+q[Nx-1,0])*(u_1[Nx,0]-u_1[Nx-1,0]))*(1./dy**2)
    delta = Dxx+Dyy
    u[Nx,0] = ( (2*k1*dt)/(k1-1)  )*V_[Nx,0]+r2*u_1[Nx,0]+r3*(delta +f_[Nx,0])
    
    
    #i = Nx and 0<= j < Ny
    Dxx = (0.5*(q[Nx,0:Ny]+q[Nx-1,0:Ny])*(u_1[Nx-1,0:Ny]-u_1[Nx,0:Ny])-0.5*(q[Nx,0:Ny]+q[Nx-1,0:Ny])*(u_1[Nx,0:Ny]-u_1[Nx-1,0:Ny]))*(1./dx**2)
    Dyy = (0.5*(q[Nx,0:Ny]+q[Nx,1:Ny+1])*(u_1[Nx,1:Ny+1]-u_1[Nx,0:Ny])-0.5*(q[Nx,0:Ny]+q[Nx,1:Ny+1])*(u_1[Nx,0:Ny]-u_1[Nx,1:Ny+1]))*(1./dy**2)
    delta = Dxx+Dyy
    u[Nx,0:Ny] = ( (2*k1*dt)/(k1-1)  )*V_[Nx,0:Ny]+r2*u_1[Nx,0:Ny]+r3*(delta + f_[Nx,0:Ny])
    # i = Nx, j = Ny
    Dxx = (0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx-1,Ny]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx,Ny]-u_1[Nx-1,Ny]))*(1./dx**2)
    Dyy = (0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny-1]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny]-u_1[Nx,Ny-1]))*(1./dy**2)
    delta = Dxx+Dyy
    u[Nx,Ny] = ( (2*k1*dt)/(k1-1)  )*V_[Nx,Ny]+r2*u_1[Nx,Ny]+r3*(delta + f_[Nx,Ny] )
    
    
    #j = Nx and 0<= i < Nx
    Dxx = (0.5*(q[0:Nx,Ny]+q[0:Nx,Ny-1])*(u_1[0:Nx,Ny-1]-u_1[0:Nx,Ny])-0.5*(q[0:Nx,Ny]+q[0:Nx,Ny-1])*(u_1[0:Nx,Ny]-u_1[0:Nx,Ny-1]))*(1./dx**2)
    Dyy = (0.5*(q[0:Nx,Ny]+q[1:Nx+1,Ny])*(u_1[1:Nx+1,Ny]-u_1[0:Nx,Ny])-0.5*(q[0:Nx,Ny]+q[1:Nx+1,Ny])*(u_1[0:Nx,Ny]-u_1[1:Nx+1,Ny]))*(1./dy**2)
    delta = Dxx+Dyy
    u[0:Nx,Ny] = ( (2*k1*dt)/(k1-1)  )*V_[0:Nx,Ny]+r2*u_1[0:Nx,Ny]+r3*(delta + f_[0:Nx,Ny])
    # j = Ny, i = Nx
    Dxx = (0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny-1]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny]-u_1[Nx,Ny-1]))*(1./dx**2)
    Dyy = (0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx-1,Ny]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx,Ny]-u_1[Nx-1,Ny]))*(1./dy**2)
    delta = Dxx+Dyy
    u[Nx,Ny] = ( (2*k1*dt)/(k1-1)  )*V_[Nx,Ny]+r2*u_1[Nx,Ny]+r3*(delta + f_[Nx,Ny])
        
    
    # switch variable before next step
   
    u_2[:] = u_1 ; u_1[:] = u
    
    # solution at subsequent step
    Nt = int(round(T/dt))
    for n in range(Nt):
        rx,ry = (1./dx**2),(1./dy**2)
        Dxx = (0.5*(q[1:Nx,1:Ny]+q[2:Nx+1,1:Ny])*(u_1[2:Nx+1,1:Ny]-u_1[1:Nx,1:Ny])-0.5*(q[1:Nx,1:Ny]+q[0:Nx-1,1:Ny])*(u_1[1:Nx,1:Ny]-u_1[0:Nx-1,1:Ny]))*rx
        Dyy = (0.5*(q[1:Nx,1:Ny]+q[1:Nx,2:Ny+1])*(u_1[1:Nx,2:Ny+1]-u_1[1:Nx,1:Ny])-0.5*(q[1:Nx,1:Ny]+q[1:Nx,0:Ny-1])*(u_1[1:Nx,1:Ny]-u_1[1:Nx,0:Ny-1]))*ry
        delta = Dxx+Dyy
        u[1:Nx,1:Ny] = k1*u_2[1:Nx,1:Ny]+k2*u_1[1:Nx,1:Ny]+k3*(delta + f_[1:Nx,1:Ny])
                
        # boundary condition

        # i = 0, 0<=j<Ny
        Dxx = (0.5*(q[0,0:Ny]+q[1,0:Ny])*(u_1[1,0:Ny]-u_1[0,0:Ny])-0.5*(q[0,0:Ny]+q[1,0:Ny])*(u_1[0,0:Ny]-u_1[1,0:Ny]))*(1./dx**2)
        Dyy = (0.5*(q[0,0:Ny]+q[0,1:Ny+1])*(u_1[0,1:Ny+1]-u_1[0,0:Ny])-0.5*(q[0,0:Ny]+q[0,1:Ny+1])*(u_1[0,0:Ny]-u_1[0,1:Ny+1]))*(1./dy**2)
        delta = Dxx+Dyy
        u[0,0:Ny] = k1*u_2[0,0:Ny]+k2*u_1[0,Ny]+k3*(delta + f_[0,0:Ny])
        # i = 0 and  j = Ny
        Dxx = (0.5*(q[0,Ny]+q[1,Ny])*(u_1[1,Ny]-u_1[0,Ny])-0.5*(q[0,Ny]+q[1,Ny])*(u_1[0,Ny]-u_1[1,Ny]))*(1./dx**2)
        Dyy = (0.5*(q[0,Ny]+q[0,Ny-1])*(u_1[0,Ny-1]-u_1[0,Ny])-0.5*(q[0,Ny]+q[0,Ny-1])*(u_1[0,Ny]-u_1[0,Ny-1]))*(1./dy**2)
        delta = Dxx+Dyy
        u[0,Ny] = k1*u_2[0,Ny]+k2*u_1[0,Ny]+k3*(delta + f_[0,Ny]) 
        
                
        #j = 0 and 0<= i < Ny
        Dxx = (0.5*(q[0:Nx,0]+q[0:Nx,1])*(u_1[0:Nx,1]-u_1[0:Nx,0])-0.5*(q[0:Nx,0]+q[0:Nx,1])*(u_1[0:Nx,0]-u_1[0:Nx,1]))*(1./dx**2)
        Dyy = (0.5*(q[0:Nx,0]+q[1:Nx+1,0])*(u_1[1:Nx+1,0]-u_1[0:Nx,0])-0.5*(q[0:Nx,0]+q[1:Nx+1,0])*(u_1[0:Nx,0]-u_1[1:Nx+1,0]))*(1./dy**2)
        delta = Dxx+Dyy
        u[0:Nx,0] = k1*u_2[0:Nx,0]+k2*u_1[0:Nx,0]+k3*(delta + f_[0:Nx,0])
        # j = 0 and  i = Nx
        Dxx = (0.5*(q[Nx,0]+q[Nx,1])*(u_1[Nx,1]-u_1[Nx,0])-0.5*(q[Nx,0]+q[Nx,1])*(u_1[Nx,0]-u_1[Nx,1]))*(1./dx**2)
        Dyy = (0.5*(q[Nx,0]+q[Nx-1,0])*(u_1[Nx-1,0]-u_1[Nx,0])-0.5*(q[Nx,0]+q[Nx-1,0])*(u_1[Nx,0]-u_1[Nx-1,0]))*(1./dy**2)
        delta = Dxx+Dyy
        u[Nx,0] = k1*u_2[Nx,0]+k2*u_1[Nx,0]+k3*(delta + f_[Nx,0])
        
     
        #i = Nx and 0<= j < Ny
        Dxx = (0.5*(q[Nx,0:Ny]+q[Nx-1,0:Ny])*(u_1[Nx-1,0:Ny]-u_1[Nx,0:Ny])-0.5*(q[Nx,0:Ny]+q[Nx-1,0:Ny])*(u_1[Nx,0:Ny]-u_1[Nx-1,0:Ny]))*(1./dx**2)
        Dyy = (0.5*(q[Nx,0:Ny]+q[Nx,1:Ny+1])*(u_1[Nx,1:Ny+1]-u_1[Nx,0:Ny])-0.5*(q[Nx,0:Ny]+q[Nx,1:Ny+1])*(u_1[Nx,0:Ny]-u_1[Nx,1:Ny+1]))*(1./dy**2)
        delta = Dxx+Dyy
        u[Nx,0:Ny] = k1*u_2[Nx,0:Ny]+k2*u_1[Nx,0:Ny]+k3*(delta + f_[Nx,0:Ny])
        # i = Nx, j = Ny
        Dxx = (0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx-1,Ny]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx,Ny]-u_1[Nx-1,Ny]))*(1./dx**2)
        Dyy = (0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny-1]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny]-u_1[Nx,Ny-1]))*(1./dy**2)
        delta = Dxx+Dyy
        u[Nx,Ny] = k1*u_2[Nx,Ny]+k2*u_1[Nx,Ny]+k3*(delta + f_[Nx,Ny])
                   

        #j = Ny and 0<= i < Nx
        Dxx = (0.5*(q[0:Nx,Ny]+q[0:Nx,Ny-1])*(u_1[0:Nx,Ny-1]-u_1[0:Nx,Ny])-0.5*(q[0:Nx,Ny]+q[0:Nx,Ny-1])*(u_1[0:Nx,Ny]-u_1[0:Nx,Ny-1]))*(1./dx**2)
        Dyy = (0.5*(q[0:Nx,Ny]+q[1:Nx+1,Ny])*(u_1[1:Nx+1,Ny]-u_1[0:Nx,Ny])-0.5*(q[0:Nx,Ny]+q[1:Nx+1,Ny])*(u_1[0:Nx,Ny]-u_1[1:Nx+1,Ny]))*(1./dy**2)
        delta = Dxx+Dyy
        u[0:Nx,Ny] = k1*u_2[0:Nx,Ny]+k2*u_1[0:Nx,Ny]+k3*(delta + f_[0:Nx,Ny])
        # j = Ny, i = Nx
        Dxx = (0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny-1]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx,Ny-1])*(u_1[Nx,Ny]-u_1[Nx,Ny-1]))*(1./dx**2)
        Dyy = (0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx-1,Ny]-u_1[Nx,Ny])-0.5*(q[Nx,Ny]+q[Nx-1,Ny])*(u_1[Nx,Ny]-u_1[Nx-1,Ny]))*(1./dy**2)
        delta = Dxx+Dyy
        u[Nx,Ny] =k1*u_2[Nx,Ny]+k2*u_1[Nx,Ny]+k3*(delta + f_[Nx,Ny])
        
        # update solution for next time step
        u_2[:], u_1[:], = u_1, u
        
    cpu_time = t0 - time.clock()
    return u, x, y, cpu_time

# 
# import nose.tools as nt
# def test_constant_solution():
#     """
#     Test for u = u_constant. 
#     usage: error = test_constant_solution()
#     """
#     
#     def exact_solution():
#         return 20.9
#     
#     def I(x,y):
#         return exact_solution()
#     
#     def V(x,y):
#         return 0
#     
#     def f(x,y,t):
#         return 0
#     
#     def velocity(x,y):
#         return 1.
#     
#     Lx,Ly = 5,5
#     Nx, Ny = 2, 2
#     dt, Tt, b = 0.1, 1, 1.
#     ue = exact_solution()
#     U,r,s,cputime = vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,Tt,b)
#     difference = abs(ue-U).max()
#     nt.assert_almost_equal(difference,0,places=14)
#     return difference
#     
#     
# def test_plug():
#     """
#     Test for plug wave : I = lambda:x,y return 0 if abs(x-L/2) > value else 1
#     for a specified value.
#     compare the numericale solution with the plug wave.
#     
#     usage: error = test_plug()
#     """
#     Lx,Ly = 5,5
#     Nx, Ny = 10, 10
#     dt, T, b = 0.1, 1., 1.
#      
#     def I(x,y):
#         L = 5
#         a = 2*Nx
#         I = np.zeros((Nx+1,Ny+1))
#         x,y = x[:,0], y[0,:]
#         I[L/2.-a:L/2.+a,:] = 1
#         return I
#  
#     def V(x,y):
#         return 0
#      
#     def f(x,y,t):
#         return 0
#      
#     def velocity(x,y):
#         return 1
# 
#     
#     #r, s = np.linspace(0,Lx,Nx+1),np.linspace(0,Ly,Ny+1)
#     U,r,s,cputime = vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
#     xv, yv = r[:,np.newaxis],s[np.newaxis,:]
#     ue = np.zeros((Nx+1,Ny+1))
#     ue = I(xv,yv)
#     print ue
#      
#     difference = abs(ue-U).max()
#     nt.assert_almost_equal(difference,0,places=13)
#     return difference
#     
#     
# def test_manufactured_solution():
#     """
#     find a manufacture solution of the pde.
#     given an analytical solution u(x,y,t) find f(x,y,t).
#     compare the analytical solution with the numerical solution.
#     calculate the convergence rate p.
#     the program return te error e, anf the convergence rate p:
#     
#     e,p = test_manufacture_solution()
#     """
# #   mx,my, A, B w ,c = 0.25,0.25, 0.25.0.25,pi/4, 1
#     X,Y,T = sp.symbols("X Y T ")
#     q = X**2+Y**2
#     A, B = 0.3,0.3
#     Lx,Ly = 1,1
#     kx, ky = 5*pi/Lx, pi/Ly
#     w = sqrt(kx**2+ky**2)
#     u = (A*sp.cos(w*T) + B*sp.sin(w*T))*sp.exp(-T)*sp.cos(kx*X)*sp.sin(ky*Y)
#     
#     Dxx, Dyy = sp.diff(q*sp.diff(u,X),X), sp.diff(q*sp.diff(u,Y),Y)
#     Dt, Dtt = sp.diff(u,T), sp.diff(sp.diff(u,T),T)
#     f1 = Dtt+Dt-Dxx-Dyy
#        
#     f2 = sp.lambdify( (X,Y,T),f1, modules = "numpy" )
#       
#     def f(x,y,t):
#         return f2(x,y,t)
#             
#     def velocity(x,y):
#         return x**2 + y**2
#       
#     ut = sp.lambdify((X,Y,T),Dt,modules = "numpy")
#      
#     
#     def exact_solution(x,y,t):
#         return  (A*cos(w*t) + B*sp.sin(w*t))*exp(-t)*cos(kx*x)*sin(ky*y)
#          
#       
#     def V(x,y):
#         return ut(x,y,0)
#     
#     def I(x,y):
#         return exact_solution(x,y,0)
#     
#     Lx,Ly = 1,1
#     Nx, Ny = 5, 5
#     dt, T, b = 0.01, 0.1,1.
#     r,s = np.linspace(0,Lx,Nx+1),np.linspace(0,Lx,Nx+1)
#     xv, yv = r[:,np.newaxis],s[np.newaxis,:]
#     
#     U,r,s,cputime = vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
#     xv, yv = r[:,np.newaxis],s[np.newaxis,:]
#     ue = np.zeros((Nx+1,Ny+1))
#     ue[:,:] = I(xv,yv)
# 
#      
#     difference = abs(ue-U).max()
#     #nt.assert_almost_equal(difference,0,places=13)
#     
#     h = [0.1,0.2,0.3,0.4,0.5]
#     E = []
#     P = []
#     for dt in h:
#         U,r,s,cputime = vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
#         e = abs(U-ue).max()
#         E.append(e)
#     from math import log
# 
#     for i in range(len(E)-1):
#         p =  log( E[i]/E[i+1]) /log(h[i+1]/h[i]) 
#         P.append(abs(p))
#     print difference,P
#     

def plot():
    """
    plot solution for T = 0.1, 1,3
    """
    def I(x,y):
        return np.exp(-0.5*(x-10/2.)**2  -0.5*(y-Ly/2.)**2)
    
    def f(x,y,t):
        return 1
    
    def velocity(x,y):
        return 1.
    
    def V(x,y):
        return 1.
    Lx,Ly = 10,10
    Nx, Ny = 100, 100
    dt, T, b = 0.001, 4., 1.
    r,s = np.linspace(0,Lx,Nx+1),np.linspace(0,Lx,Nx+1)
    for T in [4,2,1,0.5]:
        U,r,s,cputime = vectorized_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
        X, Y = meshgrid(r,s)
        fig = pl.figure(figsize=(20,20))
        ax = fig.add_subplot(1,1,1,projection='3d')
        p = ax.plot_surface(X,Y,U, rstride =1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False )
        ax.set_zlim3d(0, 3*pi);
    
        cb = fig.colorbar(p)
    pl.show()

   
if __name__=="__main__":
    plot()
    
    
    


