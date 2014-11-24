from dolfin import *
import numpy as np
from math import log
import pylab as pl
def diffusion(dt,size):
    
    # Create mesh and define function space
    #mesh = UnitSquareMesh(size, size)
    mesh = UnitIntervalMesh(size)
    V = FunctionSpace(mesh, "Lagrange", 1)
    rho = 1
    dt = 0.01
    t = 0
    f = Expression('rho*pow(x[0],2)*(-2*x[0] + 3)/6 - \
        (-12*t*x[0] + 3*t*(-2*x[0] + 3))*(pow(x[0],4)*pow((-dt + t),2)*pow((-2*x[0] + 3),2) + 36)/324 - \
        (-6*t*pow(x[0],2) + 6*t*x[0]*(-2*x[0] + 3))*(36*pow(x[0],4)*pow((-dt + t),2)*(2*x[0] - 3) + \
        36*pow(x[0],3)*pow((-dt + t),2)*pow((-2*x[0] + 3),2))/5832',
        rho=rho, t=t, dt=dt)

    
    u = TrialFunction(V)
    v = TestFunction(V)
    
    g = Constant(0.0)
    T = 0.1
    rho = 1
    I = Expression('t*x[0]*x[0]*(0.5-(1./3)*x[0])',t=0)
    
    def alpha(u):
        return 1+u**2
    u_1 = interpolate(I,V) 
    #define trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)
        
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
    
    #plot(u)
    u_a = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)', t=0.1)
    
    u_e = interpolate(u_a, V)
    Ue = u_e.vector().array()
    
    U = u.vector().array()
    x = np.linspace(0,1,len(U))
    
    pl.plot(x,U,'-',x,Ue,'ro')
    pl.show()
    #pl.ylim([0,0.008])
    #pl.xlim([0,1])
    #plot(u)
    #plot(u,interactive=True)
    #plot(u_e,interactive=True)
    
    e = np.abs(u_e.vector().array()-u.vector().array()).max()
      
    L_2norm = np.sqrt(np.sum(e**2)/u.vector().array().size)
    return L_2norm
    

 
  
dt = 0.01
s = 2000
print diffusion(dt,s)



    
    





