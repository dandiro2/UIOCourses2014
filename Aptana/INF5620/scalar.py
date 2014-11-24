import numpy as np
from math import sqrt, sin, exp
from pylab import *
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D
from scitools.std import *

def wave(Nx,Ny,Lx,Ly,V,b,f,I,T,s):
    
    q = 10
    dx, dy = float(Lx/Nx), float(Ly/Ny)
    
    rx, ry = (1./dx)**2, (1./dy)**2
    
    dt = (s/q)*sqrt(1./(rx+ry))
    #dt = q*dx
    
    k1, k2, k3 = (dt*b*0.5-1)/(1+dt*b*0.5) , 2./(1+dt*b*0.5), (dt**2)/(1+dt*b*0.5)
    
    u  = np.zeros((Nx+1,Ny+1)) # solution array at the new time level
    u_1 = np.zeros((Nx+1,Ny+1))# solution at 1 time level back
    u_2 = np.zeros((Nx+1,Ny+1))# solution at 2 time level back
    
    x, y  = np.linspace(-Lx,Lx,Nx+1), np.linspace(-Ly,Ly,Ny+1)
    

    
    # load initial condition into u_1
    #def I(X,Y):
        #return exp(-0.5*(X-Lx/2.)**2  - 0.5*(Y-Ly/2.)**2)
        #return 2*sin(pi*0.25*X)*sin(pi*0.25*Y)
        #return 20.
    
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            u_1[i,j] = I(x[i],y[j])
            
#     # special formula for first time step
    for i in range(1,Nx):
        for j in range(1,Ny):
            Dxx = (q/dx**2)*(   u_1[i+1,j]-2*u_1[i,j] + u_1[i-1,j])
            Dyy = (q/dy**2)*(   u_1[i,j+1]-2*u_1[i,j] + u_1[i,j-1])
            delt = Dxx+Dyy
              
            u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delt+f(x[i],y[j])) )
              
#     # inssert boundary condition
    i = 0
    for j in range(0,Ny+1):
        if j == Ny:
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
        else:
            Dyy =(2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
   
        Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f(x[i],y[j])) )
           
#          
    j = 0
    for i in range(0,Nx+1):
        if i == Nx:
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
        else:
            Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
  
        Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f(x[i],y[j])) )
          
    i = Nx
    for j in range(0,Ny+1):
          
        if j == 0:
            Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
        else:
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
  
        Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f(x[i],y[j])) )
          
          
    j = Ny
    for i in range(0,Nx+1):
          
        if i == 0:
            Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
        else:
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
  
        Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
        delta0 = Dxx+Dyy
        u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f(x[i],y[j])) )
         
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
                u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta+f(x[i],y[j]))

         
        
                 
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
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f(x[i],y[j]))
  
         
         
        
        for i in range(0,Nx+1):
            j = 0
            if i == Nx:
                Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
            else:
                Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
  
            Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_2[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f(x[i],y[j]))
  
          
        
        for j in range(0,Ny+1):
            i = Nx
            if j == 0:
                Dyy = (2*q/dx**2)*(u_1[i,j+1]-u_1[i,j])
            else:
                Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
              
            Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_1[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f(x[i],y[j]))
  
          
          
    
        for i in range(0,Nx+1):
            j = Ny
            if i == 0:
                Dxx = (2*q/dx**2)*(u_1[i+1,j]-u_1[i,j])
            else:
                Dxx = (2*q/dx**2)*(u_1[i-1,j]-u_1[i,j])
  
            Dyy = (2*q/dx**2)*(u_1[i,j-1]-u_1[i,j])
            delta0 = Dxx+Dyy
            #u[i,j] = (1/(1-k1))*( k2*u_2[i,j]-2*k1*dt*V +k3*(delta0+f) )
            u[i,j] = k1*u_2[i,j]+k2*u_1[i,j]+k3*(delta0+f(x[i],y[j]))
  
        u_2[:], u_1[:], = u_1, u
                 

    return u


import nose.tools as nt
def test_constant_solution():
    """
    Test for u = u_constant.
    """
    
    def exact_solution():
        return 7900
    
    def I(x,y):
        return exact_solution()
    
    def f(x,y):
        return 0
    
    S = [0.81,0.9,1]
    for s in S:    
        Nx,Ny,Lx,Ly,V,b,T,s = 100,100,5,5,0,1.,0.5,s
    
        ue = exact_solution()
        U,r,e = wave(Nx,Ny,Lx,Ly,V,b,f,I,T,s)
        difference = abs(ue-U).max()

        #nt.assert_almost_equal(difference,0,places=14)
        tol = 1e-15
        
        if difference <= tol:
            break
    print difference
    print "s= ",s
    nt.assert_almost_equal(difference,0,places=14)

def test_plug_wave_solution():
    """
    Test for u = u_constant.
    """
    s = 0.5
    
    Nx,Ny,Lx,Ly,V,b,T,s = 100,100,5.,5.,0.,1.,0.5,s
    
    exact_solutionx = lambda x,y: 0 if abs( (x-Lx/2.0) ) > 0.1 else 1
    
    I = exact_solutionx
    
    def f(x,y):
        return 0

    U = wave(Nx,Ny,Lx,Ly,V,b,f,I,T,s)
    

     
    x,y = np.linspace(-Lx,Lx,Nx+1),np.linspace(-Ly,Lx,Ny+1)
    ue = np.zeros((Nx+1,Ny+1))
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            ue = exact_solutionx(x[i],y[j])
                
    difference = abs(ue-U).max()
 
        
    print difference
    print "s= ",s
    nt.assert_almost_equal(difference,0,places=14)

   
    
    
    
    
    
    
    
def plot():
    U, x, y = wave()
    X, Y = meshgrid(x,y)
    
    fig = pl.figure(figsize=(14,16))
    #ax = fig.add_subplot(1,1,1,projection='3d')
    #p = ax.plot_surface(X,Y,U, rstride =1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False )
    
    ax = fig.add_subplot(1,1,1, projection='3d')
    p = ax.plot_surface(X,Y,U, rstride =5, cstride=5, cmap=cm.coolwarm, linewidth=0, antialiased=False )

    ax.plot_surface(X, Y, U, rstride=5, cstride=5, alpha=0.25)
    cset = ax.contour(X, Y, U, zdir='z', offset=-pi, cmap=cm.coolwarm)
    #cset = ax.contour(X, Y, U, zdir='x', offset=-pi, cmap=cm.coolwarm)
    #cset = ax.contour(X, Y, U, zdir='y', offset=3*pi, cmap=cm.coolwarm)
    #
    #ax.set_xlim3d(0, 2*pi);
    #ax.set_ylim3d(0, 3*pi);
    ax.set_zlim3d(0, 3*pi);
    
    cb = fig.colorbar(p)
    pl.show()
    
    
    
#plot()
test_plug_wave_solution()   




