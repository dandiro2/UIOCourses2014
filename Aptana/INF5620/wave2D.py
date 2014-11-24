import numpy as np
from math import sqrt, sin, exp
from pylab import *
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

def wave():
    
    Nx, Ny = 100, 100
    Lx, Ly = 5.,5.
    s = 0.9
    q = 1.
    b = 1.
    V = 0
    f = 0
    T = 1
    dx, dy = float(Lx/Nx), float(Ly/Ny)
    
    rx, ry = (1./dx)**2, (1./dy)**2
    
    dt = (s/q)*sqrt(1./(rx+ry))
    
    k1, k2, k3 = (dt*b*0.5-1)/(1+dt*b*0.5) , 2./(1+dt*b*0.5), (dt**2)/(1+dt*b*0.5)
    
    u  = np.zeros((Nx+1,Ny+1)) # solution array at the new time level
    u_1 = np.zeros((Nx+1,Ny+1))# solution at 1 time level back
    u_2 = np.zeros((Nx+1,Ny+1))# solution at 2 time level back
    
    x, y  = np.linspace(-Lx,Lx,Nx+1), np.linspace(-Ly,Ly,Ny+1)
    

    
    # load initial condition into u_1
    def I(X,Y):
        #return exp(-0.5*(X-Lx/2.)**2  - 0.5*(Y-Ly/2.)**2)
        #return 2*sin(pi*0.25*X)*sin(pi*0.25*Y)
        return 200

    for i in range(0,Nx+1):
        for j in range(Ny+1):
            u_1[i,j] = I(x[i],y[j])
            
#     # special formula for first time step
    for i in range(1,Nx):
        for j in range(1,Ny):
            Dxx = (q/dx**2)*(   u_1[i+1,j]-2*u_1[i,j] + u_1[i-1,j])
            Dyy = (q/dy**2)*(   u_1[i,j+1]-2*u_1[i,j] + u_1[i,j-1])
            delt = Dxx+Dyy
              
            u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delt+f) )
              
#     # inssert boundary condition
    i = 0
    for j in range(0,Ny+1):
        if j == Ny:
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
        else:
            Dyy =(2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
   
        Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
           
#          
    j = 0
    for i in range(0,Nx+1):
        if i == Nx:
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
        else:
            Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
  
        Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
          
    i = Nx
    for j in range(0,Ny+1):
          
        if j == 0:
            Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
        else:
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
  
        Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
          
          
    j = Ny
    for i in range(0,Nx+1):
          
        if i == 0:
            Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
        else:
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
  
        Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
         
    # switch variable before next step
    u_2[:], u_1[:] = u_1, u
    #solution at the subsequent step
    
    Nt = int(round(T/dt))
    for n in range(Nt):
      
        for i in range(1,Nx):
            for j in range(1,Ny):
                 
                Dxx = (q/dx**2)*(   u_1[i+1,j]-2*u_1[i,j] + u_1[i-1,j])
                Dyy = (q/dy**2)*(   u_1[i,j+1]-2*u_1[i,j] + u_1[i,j-1])
                delta = Dxx+Dyy
                u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta+f)

         
        
                 
        #boundary condition
        
        for j in range(0,Ny+1):
            i = 0
            if j == Ny:
                Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
            else:
                Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
   
            Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_2[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f)
  
         
         
        
        for i in range(0,Nx+1):
            j = 0
            if i == Nx:
                Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
            else:
                Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
  
            Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_2[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f)
  
          
        
        for j in range(0,Ny+1):
            i = Nx
            if j == 0:
                Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
            else:
                Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
              
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f)
  
          
          
    
        for i in range(0,Nx+1):
            j = Ny
            if i == 0:
                Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
            else:
                Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
  
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_2[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f)
  
        u_2[:], u_1[:], = u_1, u
                 

    return u,x,y


import nose.tools as nt
def test_constant_solution():
    """
    Test for u = u_constant.
    """
    
    def exact_solution():
        return 200
    
    ue = exact_solution()
    U,r,s = wave()
    difference = abs(ue-U).max()

    nt.assert_almost_equal(difference,0,places=14)
    print difference
    #tol = 1e-15
    #assert difference <= tol
    
    
def plot():

    X, Y = meshgrid(x,y)
    U = ...
    fig = pl.figure(figsize=(14,16))
    ax = fig.add_subplot(1,1,1, projection='3d')
    p = ax.plot_surface(X,Y,U, rstride =5, cstride=5, cmap=cm.coolwarm, linewidth=0, antialiased=False )
    ax.plot_surface(X, Y, U, rstride=5, cstride=5, alpha=0.25)
    cset = ax.contour(X, Y, U, zdir='z', offset=-pi, cmap=cm.coolwarm)
    ax.set_zlim3d(0, 3*pi);
    cb = fig.colorbar(p)
    pl.show()
    
    
    
#plot()
test_constant_solution()   




