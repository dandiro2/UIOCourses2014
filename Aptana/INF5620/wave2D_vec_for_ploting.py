import numpy as np
from math import sqrt, sin, exp
from pylab import *
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

def wave():
    
    Nx, Ny = 200, 200
    Lx, Ly = 10.,10.
    s = 0.9
    q = 1.
    b = 1
    V = 1.
    f = 0.1
    T = 4
    start = 0
    
    dx, dy = float(Lx/Nx), float(Ly/Ny)
    
    rx, ry = (1./dx)**2, (1./dy)**2
    
    dt = (s/q)*sqrt(1./(rx+ry))
    
    k1, k2, k3 = (dt*b*0.5-1)/(1+dt*b*0.5) , 2./(1+dt*b*0.5), (dt**2)/(1+dt*b*0.5)
    
    u  = np.zeros((Nx+1,Ny+1)) # solution array at the new time level
    u_1 = np.zeros((Nx+1,Ny+1))# solution at 1 time level back
    u_2 = np.zeros((Nx+1,Ny+1))# solution at 2 time level back
    
    x, y  = np.linspace(start,Lx,Nx+1), np.linspace(start,Ly,Ny+1)
    
    xv, yv = x[:,newaxis],y[newaxis,:]
    
    #load initial condition into u_1
    def I(X,Y):
        return np.exp(-0.5*(X-Lx/2.)**2  -0.5*(Y-Ly/2.)**2)
        #return 2*sin(pi*0.25*X)*sin(pi*0.25*Y)
        #return 20.

#     for i in range(0,Nx+1):
#         for j in range(0,Ny+1):
#             u[i,j] = I(x[i],y[j])
            
    
    u_1[:,:] = I(xv,yv)
#             
# #     # special formula for first time step
# 
    Dxx = (q/dx**2)*(   u_1[2:Nx+1,1:Ny]-2*u_1[1:Nx,1:Nx] + u_1[0:Nx-1,1:Ny])
    Dyy = (q/dy**2)*(   u_1[1:Nx,2:Nx+1]-2*u_1[1:Nx,1:Ny] + u_1[1:Nx,0:Ny-1])
    delt = Dxx+Dyy      
    u[1:Nx,1:Ny] = (1/(1-k1))*( k2*u_1[1:Nx,1:Ny]-2*k1*dt*V +k3*(delt+f) )
               
    # inssert boundary condition
 
    #######i = 0; 0<= j < Ny#############
    Dxx = (2*q/dx**2)*(u_1[1,0:Ny]-u_1[0,0:Ny])
    Dyy =(2*q/dx**2)*(u_1[0,1:Ny+1]-u_1[0,0:Ny])
    delta0 = Dxx+Dyy
    u[0,0:Ny] = (1/(1-k1))*( k2*u_1[0,0:Ny]-2*k1*dt*V +k3*(delta0+f) )          
     
    # i = 0; j = Ny
    Dyy = (2*q/dx**2)*(u_1[0,Ny-1]-u_1[0,Ny])
    Dxx = (2*q/dx**2)*(u_1[1,Ny]-u_1[0,Ny])
    delta0 = Dxx+Dyy
    u[0,Ny] = (1/(1-k1))*( k2*u_1[0,Ny]-2*k1*dt*V +k3*(delta0+f) )
     
 
    #########j = 0; 0<= i < Nx###########
    Dxx = (2*q/dx**2)*(u_1[0:Nx,1]-u_1[0:Ny,0])
    Dyy =(2*q/dx**2)*(u_1[1:Nx+1,0]-u_1[0:Nx,0])
    delta0 = Dxx+Dyy
    u[0:Nx,0] = (1/(1-k1))*( k2*u_1[0:Nx,0]-2*k1*dt*V +k3*(delta0+f) )          
     
    # j = 0; i = Nx
    Dyy = (2*q/dx**2)*(u_1[Nx-1,0]-u_1[Nx,0])
    Dxx = (2*q/dx**2)*(u_1[Nx,1]-u_1[Nx,0])
    delta0 = Dxx+Dyy
    u[Nx,0] = (1/(1-k1))*( k2*u_1[Nx,0]-2*k1*dt*V +k3*(delta0+f) )
      
      
      
 
    ######i = Nx, 1 <= j <Ny+1########
    Dxx = (2*q/dx**2)*(u_1[Nx-1,1:Ny+1]-u_1[Nx,1:Ny+1])
    Dyy = (2*q/dx**2)*(u_1[Nx,0:Ny]-u_1[Nx,1:Ny+1])
    delta0 = Dxx+Dyy
    u[Nx,1:Ny+1] = (1/(1-k1))*( k2*u_1[Nx,1:Ny+1]-2*k1*dt*V +k3*(delta0+f) )
        
    # i = Nx; j = 0   
    Dxx = (2*q/dx**2)*(u_1[Nx-1,0]-u_1[Nx,0])
    Dyy = (2*q/dx**2)*(u_1[Nx,1]-u_1[Nx,0])
    delta0 = Dxx+Dyy
    u[Nx,0] = (1/(1-k1))*( k2*u_1[Nx,0]-2*k1*dt*V +k3*(delta0+f) )
           
           
           
           
 
    #######j = Ny, 1 <= i <Nx+1#####
    Dxx = (2*q/dx**2)*(u_1[1:Nx+1,Ny-1]-u_1[1:Nx+1,Ny])
    Dyy = (2*q/dx**2)*(u_1[0:Nx,Ny]-u_1[1:Nx+1,Ny])
    delta0 = Dxx+Dyy
    u[1:Nx+1,Ny] = (1/(1-k1))*( k2*u_1[1:Nx+1,Ny]-2*k1*dt*V +k3*(delta0+f) )
        
    # j = Ny; i = 0   
    Dxx = (2*q/dx**2)*(u_1[0,Ny-1]-u_1[0,Ny])
    Dyy = (2*q/dx**2)*(u_1[1,Ny]-u_1[0,Ny])
    delta0 = Dxx+Dyy
    u[0,Ny] = (1/(1-k1))*( k2*u_1[0,Ny]-2*k1*dt*V +k3*(delta0+f) )
           
          
    # switch variable before next step
    u_2[:], u_1[:] = u_1, u
     
    #solution at the subsequent step
    Nt = int(round(T/dt))
    for n in range(Nt):
 
        Dxx = (q/dx**2)*(   u_1[2:Nx+1,1:Nx]-2*u_1[1:Nx,1:Nx] + u_1[0:Nx-1,1:Nx])
        Dyy = (q/dy**2)*(   u_1[1:Nx,2:Ny+1]-2*u_1[1:Nx,1:Ny] + u_1[1:Nx,0:Ny-1])
        delta = Dxx+Dyy
        u[1:Nx,1:Ny] = k1*u_2[1:Nx,1:Ny]+k2*u_1[1:Nx,1:Ny]+k3*(delta+f)
         
                  
        #boundary condition
         
        ### i = 0, 0<= j <Ny
        Dxx = (2*q/dx**2)*(u_1[1,0:Ny]-u_1[0,0:Ny])
        Dyy = (2*q/dx**2)*(u_1[0,1:Ny+1]-u_1[0,0:Ny])
        delta0 = Dxx+Dyy
        u[0,0:Ny] = k1*u_2[0,0:Ny]+k2*u_1[0,0:Ny]+k3*(delta0+f)
        # i = 0, J = Ny
        Dxx = (2*q/dx**2)*(u_1[1,Ny]-u_1[0,Ny])
        Dyy = (2*q/dx**2)*(u_1[0,Ny-1]-u_1[0,Ny])
        delta0 = Dxx+Dyy
        u[0,Ny] = k1*u_2[0,Ny]+k2*u_1[0,Ny]+k3*(delta0+f)
         
         
 
        ### j = 0, 0<= i <Nx
        Dxx = (2*q/dx**2)*(u_1[0:Nx,1]-u_1[0:Nx,0])
        Dyy = (2*q/dx**2)*(u_1[1:Nx+1,0]-u_1[0:Nx,0])
        delta0 = Dxx+Dyy
        u[0:Nx,0] = k1*u_2[0:Nx,0]+k2*u_1[0:Nx,0]+k3*(delta0+f)
        # j = 0, i = Nx
        Dxx = (2*q/dx**2)*(u_1[Nx,1]-u_1[Nx,0])
        Dyy = (2*q/dx**2)*(u_1[Nx-1,0]-u_1[Nx,0])
        delta0 = Dxx+Dyy
        u[Nx,0] = k1*u_2[Nx,0]+k2*u_1[Nx,0]+k3*(delta0+f)
         
         
         
        ##### i = Nx; 1<= j <Ny+1######
        Dyy = (2*q/dx**2)*(u_1[Nx,0:Ny]-u_1[Nx,1:Ny+1])
        Dxx = (2*q/dx**2)*(u_1[Nx-1,1:Ny+1]-u_1[Nx,1:Ny+1])
        delta0 = Dxx+Dyy
        u[Nx,1:Ny+1] = k1*u_2[Nx,1:Ny+1]+k2*u_1[Nx,1:Ny+1]+k3*(delta0+f)
        #i = Nx; j = 0
        Dxx = (2*q/dx**2)*(u_1[Nx-1,0]-u_1[Nx,0])
        Dyy = (2*q/dx**2)*(u_1[Nx,1]-u_1[Nx,0])
        delta0 = Dxx+Dyy
        u[Nx,0] = k1*u_2[Nx,0]+k2*u_1[Nx,0]+k3*(delta0+f)
           
 
 
        ##### j = Ny; 1<= i <Nx+1######
        Dyy = (2*q/dx**2)*(u_1[0:Nx,Ny]-u_1[1:Nx+1,Ny])
        Dxx = (2*q/dx**2)*(u_1[1:Nx+1,Ny-1]-u_1[1:Nx+1,Ny])
        delta0 = Dxx+Dyy
        u[1:Nx+1,Ny] = k1*u_2[1:Nx+1,Ny]+k2*u_1[1:Nx+1,Ny]+k3*(delta0+f)
        #j = Ny; i = 0
        Dxx = (2*q/dx**2)*(u_1[0,Ny-1]-u_1[0,Nx])
        Dyy = (2*q/dx**2)*(u_1[1,Ny]-u_1[0,Ny])
        delta0 = Dxx+Dyy
        u[0,Ny] = k1*u_2[0,Ny]+k2*u_1[0,Ny]+k3*(delta0+f)
         
        # update for next time step
        u_2[:], u_1[:], = u_1, u
                  

    return u,x,y
    
    


import nose.tools as nt
def test_constant_solution():
    """
    Test for u = u_constant.
    """
    
    def exact_solution():
        return 20.
    
    ue = exact_solution()
    U,r,s = wave()
    difference = abs(ue-U).max()

    nt.assert_almost_equal(difference,0,places=14)
    print difference
    #tol = 1e-15
    #assert difference <= tol
    
    
def plot():
    U, x, y = wave()
    X, Y = meshgrid(x,y)
    fig = pl.figure(figsize=(20,20))
    ax = fig.add_subplot(1,1,1,projection='3d')
    p = ax.plot_surface(X,Y,U, rstride =1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False )
    ax.set_zlim3d(0, 3*pi);
    
    cb = fig.colorbar(p)
    pl.show()

    
    
    
plot()





