import numpy as np
import matplotlib
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scitools.std import movie,glob
import os

# remove old files of given format
for filename in glob('wave*.png'):
    os.remove(filename)


    
def solve(I,q,f,V,b,dt,dx,dy,Lx,Ly,T):
    
    #dt=0.02;dx=0.05;dy=0.05;T=1.0;Lx=1.0;Ly=1.0
    TT = int(round(T/dt))+1
    Nx = int(round(Lx/dx))+1
    Ny = int(round(Ly/dy))+1
    t = np.linspace(0,T,TT)
    x = np.linspace(0,Lx,Nx)
    y = np.linspace(0,Ly,Ny)
    u = np.zeros((TT,Nx+2,Ny+2)) 
    un = np.zeros((TT,Nx+2,Ny+2)) 
    xv, yv = x[:,np.newaxis], y[np.newaxis,:]
    #X = X.transpose()transpose
    #Y = Y.transpose()
    u[0,1:-1,1:-1] = I(X,Y)     # Initial Condition

    # Boundaries
    u[0,0,1:-1] = u[0,2,1:-1]
    u[0,Nx+1,1:-1] = u[0,Nx-1,1:-1]
    u[0,1:-1,0] = u[0,1:-1,2]
    u[0,1:-1,Ny+1] = u[0,1:-1,Ny-1]


    
    rx = (dt/dx)**2
    ry = (dt/dy)**2

    dqu_dx = rx*(q(xv+dx/2,yv)*(u[0,2:,1:-1]-u[0,1:-1,1:-1])-q(xv-dx/2,yv)*(u[0,1:-1,1:-1]-u[0,0:-2,1:-1]))

    dqu_dy = ry*(q(xv,yv+dy/2)*(u[0,1:-1,2:]-u[0,1:-1,1:-1])-q(xv,yv-dy/2)*(u[0,1:-1,1:-1]-u[0,1:-1,0:-2]))


    # first step:
    u[1,1:-1,1:-1] = 0.5*(dqu_dx+dqu_dy 
                + 2*u[0,1:-1,1:-1]
                + 2*dt*V(xv,yv)*(1-dt*b/2)
                + (dt**2)*f(xv,yv,0))
    # boundaries
    u[1,0,1:-1] = u[1,2,1:-1]
    u[1,Nx+1,1:-1] = u[1,Nx-1,1:-1]
    u[1,1:-1,0] = u[1,1:-1,2]
    u[1,1:-1,Ny+1] = u[1,1:-1,Ny-1]

    # first step done, now loop to find all u[n] for n=2,...,M-1
    for n in range(1,TT-1):
          
        # expression for (q*u_x)_x
        dqu_dx = rx*(q(xv+dx/2,yv)*(u[n,2:,1:-1]-u[n,1:-1,1:-1])
                 - q(xv-dx/2,yv)*(u[n,1:-1,1:-1]-u[n,0:-2,1:-1]))    
        # expression for (q*u_y)_y
        dqu_dy = ry*(q(xv,yv+dy/2)*(u[n,1:-1,2:]-u[n,1:-1,1:-1])
                             - q(xv,yv-dy/2)*(u[n,1:-1,1:-1]-u[n,1:-1,0:-2]))
     
        # finite difference scheme for inner points
        u[n+1,1:-1,1:-1] = (1/(1+b*dt/2))*(dqu_dx+dqu_dy
                             + 2*u[n,1:-1,1:-1]
                             + u[n-1,1:-1,1:-1]*(dt*b/2 -1)
                             + (dt**2)*f(xv,yv,t[n]))

        # Boundary conditions
        u[n+1,Nx+1,1:-1] = u[n+1,Nx-1,1:-1]

        u[n+1,1:-1,0] = u[n+1,1:-1,2]
        u[n+1,1:-1,Ny+1] = u[n+1,1:-1,Ny-1]
        
        un = u[:,1:-1,1:-1]
    return un






    

    