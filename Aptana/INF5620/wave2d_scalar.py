"""
solve the 2D wave equation with variable wave velocity q(x,y).
d2u/dt2 + b(du/dt) = d/dx(q(x,y)du/dx) + d/dy(q(x,y)du/dy) +f(x,y,t)
BC: du/dn = 0 (Newman BC)
IC: du(x,y,0)/dt = V(x,y)
    u(x,y,0)     = I(x,y)
"""

import numpy as np
from math import sqrt, sin, exp,cos, pi
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D
import sympy as sp
from sympy import *

def scalar_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b):
    """"
    implementation of the scalar solver.
    I(x,y) is the initial solution, V(x,y) is the initial velocity (at t=0), q(x,y) is the 
    wave velocity at subsequent time. Lx and Ly are the length of the domain along the x and y
    axis. Nx and Ny are number of points along the x and y axis. dt is the time step, T is the 
    total simulation time. ste1BC and next_stepBC are the boundary condition at the first time
    step and the subsequent time step respectively.
    """
    
    u   = np.zeros((Nx+1,Ny+1)) # solution array
    u_1 = np.zeros((Nx+1,Ny+1)) # solution u[n]
    u_2 = np.zeros((Nx+1,Ny+1)) # solution u[n-2]
    x  = np.linspace(0, Lx, Nx+1) # mesh points in x dir
    y  = np.linspace(0, Ly, Ny+1) # mesh points in y dir
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    
    
#     stability_limit = (1/float(max(q)))*(1/sqrt(1/dx**2 + 1/dy**2))
#     if dt <= 0: # 
#         safety_factor = -dt # use negative dt as safety factor
#         dt = safety_factor*stability_limit
#     elif dt > stability_limit:
#         print 'error: dt=%g exceeds the stability limit %g' %(dt, stability_limit)
#         
    Nt = int(round(T/float(dt)))
    t = np.linspace(0, Nt*dt, Nt+1) # mesh points in time
    
    # simplification parameters k1,k2,k3,r1,r2,r3:
    k1 = ( (dt*b/2.)-1 )/( (dt*b/2.)+1 )
    k2 = 2./(  (dt*b/2.)+1   )
    k3 = (dt**2.)/((dt*b/2.)+1)
    
    
    q  = np.ones((Nx+1,Ny+1))
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            q[i,j] = velocity(x[i],y[j])
    

    r2 = k2/(1-k1)
    r3 = k3/(1-k1)
    
    import time;  t0 = time.clock()  # for measuring CPU time
    
    # initial condition
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            u_1[i,j] = I( x[i],y[j] )
    
#     # special formula for the first time step
    n = 0
    for i in range(1,Nx):
        for j in range(1,Ny):
            Dxx = ( 0.5*(q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-0.5*(q[i,j]+q[i-1,j])*(u_1[i,j]-u_1[i-1,j]))*(1./dx**2)
            Dyy = (0.5*(q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j])-0.5*(q[i,j]+q[i,j-1])*(u_1[i,j]-u_1[i,j-1]))*(1./dy**2)
            delta = Dxx+Dyy
            u[i,j] = ( (2*k1*dt)/(k1-1)  )*V( x[i],y[j] )+r2*u_1[i,j]+r3*(delta + f(x[i],y[j],t[n] ))
             
 
             
    # boundary condition
    i = 0
    for j in range(0,Ny+1):
        Ip,Jp = i+1,j+1
        Im,Jm = i-1,j-1
         
        if j == 0:
            Jm = Jp
        elif j == Ny:
            Jp = Jm
              
        Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
        Im = Ip
        Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
        delta = Dxx+Dyy
        u[i,j] = ( (2*k1*dt)/(k1-1)  )*V( x[i],y[j] )+r2*u_1[i,j]+r3*(delta + f(x[i],y[j],t[n] ))
             
    j = 0
    for i in range(0,Nx+1):
        Ip,Jp = i+1,j+1
        Im,Jm = i-1,j-1
        if i == 0:
            Im = Ip
        elif i == Nx:
            Ip = Im
        #Jm = Jp
        Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
        Jm = Jp
        Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
        delta = Dxx+Dyy
        u[i,j] = ( (2*k1*dt)/(k1-1)  )*V( x[i],y[j] )+r2*u_1[i,j]+r3*(delta + f(x[i],y[j],t[n] ))
             
    i = Nx
    for j in range(0,Ny+1):
        Ip,Jp = i+1,j+1
        Im,Jm = i-1,j-1
        if j == 0:
            Jm = Jp
        elif j == Ny:
            Jp = Jm
             
        Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
        Ip = Im
        Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
        delta = Dxx+Dyy
        u[i,j] = ( (2*k1*dt)/(k1-1)  )*V( x[i],y[j] )+r2*u_1[i,j]+r3*(delta + f(x[i],y[j],t[n] ))
             
    j = Ny
    for i in range(0,Nx+1):
        Ip,Jp = i+1,j+1
        Im,Jm = i-1,j-1
        if i == 0:
            Im = Ip
        elif i == Nx:
            Ip = Im
         
        Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
        Jp = Jm
        Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
        delta = Dxx+Dyy
        u[i,j] = ( (2*k1*dt)/(k1-1)  )*V( x[i],y[j] )+r2*u_1[i,j]+r3*(delta + f(x[i],y[j],t[n] ))
             
     
    # switch variable before next step
    #u_2[:], u_1[:] = u_1, u
    u_2[:] = u_1 ; u_1[:] = u
          
    # soution at subsequent step
    Nt = int(round(T/dt))
    for n in range(1,Nt):
       
        for i in range(1,Nx):
            for j in range(1,Ny):
                Dxx = (0.5*(q[i,j]+q[i+1,j])*(u_1[i+1,j]-u_1[i,j])-0.5*(q[i,j]+q[i-1,j])*(u_1[i,j]-u_1[i-1,j]))*(1./dx**2)
                Dyy = (0.5*(q[i,j]+q[i,j+1])*(u_1[i,j+1]-u_1[i,j])-0.5*(q[i,j]+q[i,j-1])*(u_1[i,j]-u_1[i,j-1]))*(1./dy**2)
                delta = Dxx+Dyy
                u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta + f(x[i],y[j],t[n]))
                 
        # boundary condition
        i = 0
        for j in range(0,Ny+1):
            Ip,Jp = i+1,j+1
            Im,Jm = i-1,j-1
            if j == 0:
                Jm = Jp
            elif j == Ny:
                Jp = Jm
                 
            Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
            Im = Ip
            Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
            delta = Dxx+Dyy
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta + f(x[i],y[j],t[n]))
                 
        j = 0
        for i in range(0,Nx+1):
            Ip,Jp = i+1,j+1
            Im,Jm = i-1,j-1
            if i == 0:
                Im = Ip
            elif i == Nx:
                Ip = Im
            Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
            Jm = Jp
            Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
            delta = Dxx+Dyy
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta + f(x[i],y[j],t[n]))
                 
        i = Nx
        for j in range(0,Ny+1):
            Ip,Jp = i+1,j+1
            Im,Jm = i-1,j-1
            if j == 0:
                Jm = Jp
            elif j == Ny:
                Jp = Jm
             
            Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
            Ip = Im
            Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
            delta = Dxx+Dyy
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta + f(x[i],y[j],t[n]))
                 
        j = Ny
        for i in range(0,Nx+1):
            Ip,Jp = i+1,j+1
            Im,Jm = i-1,j-1
            if i == 0:
                Im = Ip
            elif i == Nx:
                Ip = Im
             
            Dxx = (0.5*(q[i,j]+q[Ip,j])*(u_1[Ip,j]-u_1[i,j])-0.5*(q[i,j]+q[Im,j])*(u_1[i,j]-u_1[Im,j]))*(1./dx**2)
            Jp = Jm
            Dyy = (0.5*(q[i,j]+q[i,Jp])*(u_1[i,Jp]-u_1[i,j])-0.5*(q[i,j]+q[i,Jm])*(u_1[i,j]-u_1[i,Jm]))*(1./dy**2)
            delta = Dxx+Dyy
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta + f(x[i],y[j],t[n]))
                 
         
        # update solution for next time step
        u_2[:] = u_1 ; u_1[:] = u
    
    
#         
    cpu_time = t0 - time.clock()
    return u, x, y, cpu_time


import nose.tools as nt
def test_constant_solution():
    """
    Test for u = u_constant.
    """
    
    def exact_solution():
        return 21.3
    
    def I(x,y):
        return exact_solution()
    
    def V(x,y):
        return 0
    
    def f(x,y,t):
        return 0
    def velocity(x,y):
        return 1
    Lx,Ly = 5,5
    Nx, Ny = 100, 100
    dt, T, b = 0.1, 1, 1.
    ue = exact_solution()
    U,r,s,cputime = scalar_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
    difference = abs(ue-U).max()
    nt.assert_almost_equal(difference,0,places=14)
    print difference, cputime
    
def test_plug():
    """
    Test for u = u_constant.1
    """
    Lx,Ly = 5,5
    Nx, Ny = 100, 100
    dt, T, b = 0.1, 1., 1.
    
    I = lambda x,y: 0. if abs(x - (Lx/2.)) > 3 else 1.
    
    def V(x,y):
        return 0
    
    def f(x,y,t):
        return 0
    
    def velocity(x,y):
        return 1
    
    U,r,s,cputime = scalar_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
    
    ue = np.zeros((Nx+1,Ny+1))
    for i in range(Nx+1):
        for j in range(Ny+1):
            ue[i,j] = I(r[i],s[j])

    difference = abs(ue-U).max()
    nt.assert_almost_equal(difference,0,places=14)
    print difference, cputime

def test_manufactured_solution():
    """
    find a manufacture solution of the pde.
    given an analytical solution u(x,y,t) find f(x,y,t).
    find also the convergence rate
    """
    #mx,my, A, B w ,c = 0.25,0.25, 0.25.0.25,pi/4, 1
    x,y,t = symbols("x y t ")
    q = x**2+y**2
    u = (0.5*cos(pi*0.25*t) + 0.25*sin(pi*0.25*t))*exp(-t)*cos(0.25*pi*x)*sin(0.25*pi*y)

    Dxx, Dyy = diff(q*diff(u,x),x), diff(q*diff(u,y),y)
    Dt, Dtt = diff(u,t), diff(diff(u,t),t)
    f1 = Dtt+Dt-Dxx-Dyy
    
    f2 = lambdify( (x,y,t),f1, modules = "numpy" )
    
    def f(x,y,t):
        return f2(x,y,t)
        
    
    def velocity(x,y):
        return x**2+y**2
    
    ut = lambdify((x,y,t),Dt,modules = "numpy")
    def V(x,y):
        return ut(x,y,0)
    
    def exact_solution(x,y,t):
        return (0.5*cos(pi*0.25*t) + 0.25*sin(pi*0.25*t))*exp(-t)*cos(0.25*pi*x)*sin(0.25*pi*y)
    

    
    def I(x,y):
        return exact_solution(x,y,0)
    
    b = 1
    Lx,Ly, Nx, Ny = 1,1,5,5   
    dt,T = 0.01,0.1
    U,r,s,cputime = scalar_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
    
    ue = np.zeros((Nx+1,Ny+1))
    Nt = int(round(T/dt))
    
    t = np.linspace(0, Nt*dt, Nt+1) 
    
    for n in range(Nt+1):
        for i in range(Nx+1):
            for j in range(Ny+1):
                ue[i,j] = exact_solution(r[i],s[j],t[n])

    difference = abs(ue-U).max()
    
    h = [0.1,0.2,0.3,0.4,0.5]
    E = []
    for dt in h:
        U,r,s,cputime = scalar_solver(I,V,f,velocity,Lx,Ly,Nx,Ny,dt,T,b)
        e = abs(U-ue).max()
        E.append(e)
    print E
        
        
        
    #nt.assert_almost_equal(difference,0,places=14)
    print difference


    
    
#test_plug()
#test_constant_solution()
test_manufactured_solution()
