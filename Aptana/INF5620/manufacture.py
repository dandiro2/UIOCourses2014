from dolfin import *
import numpy as np
from math import log
import pylab as pl
def diffusion(dt,size):
    
    # Create mesh and define function space
    #mesh = UnitSquareMesh(size, size)
    mesh = UnitIntervalMesh(size)
    #V = FunctionSpace(mesh, "Lagrange", 1)
    #R = FunctionSpace(mesh, "R", 0)
    #W = V * R
    
    # Define variational problem
    #(u, c) = TrialFunction(W)
    #(v, d) = TestFunctions(W)
    
    f = Expression('-(x[0]*x[0]*x[0])/3. + (x[0]*x[0])/2+8*(t*t*t)*(x[0]*x[0]*x[0])/9. \
    -28*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0])/9.+7*(t*t*t)*(x[0]*x[0]*x[0]*x[0]*x[0])/2. \
    -5*(t*t*t)*(x[0]*x[0]*x[0]*x[0])/4.+ 2*t*x[0]-t',t =0)
    

    g = Constant(0.0)
    T = 0.5
    rho = 0.013
    I = Expression('t*x[0]*x[0]*(0.5-(1./3)*x[0])',t=0)

    #u_1 = interpolate(I,V) 
    
    def alpha(u):
        return 1+u**2
    
    #a = ((dt/rho)*inner(alpha(u_1)*grad(u), grad(v)) + (dt/rho)*c*v + u*d)*dx +u*v*dx
    #L = f*v*dx + g*v*ds + u_1*v*dx
    
    # Compute solution
    #w = Function(W)
    t = 0
    U = []
    while t<=T:

        V = FunctionSpace(mesh, "Lagrange", 1)
        u_1 = interpolate(I,V) 
        R = FunctionSpace(mesh, "R", 0)
        W = V * R
        w = Function(W)
    
        # Define variational problem
        (u, c) = TrialFunction(W)
        (v, d) = TestFunctions(W)
        I.t = t
        f.t = t
        a = ((dt/rho)*inner(alpha(u_1)*grad(u), grad(v)) + (dt/rho)*c*v + u*d)*dx +u*v*dx
        L = f*v*dx + g*v*ds + u_1*v*dx
    
        solve(a == L, w)
        (u, c) = w.split(True)
        u_1.assign(u)
        t+=dt
    
    plot(u)
    #u_a = Expression('t*x[0]*x[0]*(0.5-(1./3)*x[0])',t=0.05)
    
    u_e = interpolate(I, V)
    Ue = u_e.vector().array()
    
    #U = u.vector().array()
    #x = np.linspace(0,1,len(U))
    
    #pl.plot(x,U,'-',x,Ue,'r')
    #pl.show()
    #pl.ylim([0,0.008])
    #pl.xlim([0,1])
    #plot(u)
    #plot(u,interactive=True)
    plot(u_e,interactive=True)
    
    e = np.abs(u_e.vector().array()-u.vector().array()).max()
      
    L_2norm = np.sqrt(np.sum(e**2)/u.vector().array().size)
    return L_2norm
    

 
# 
# h = [0.00045,0.00040,0.00035,0.00030]
# h_hf = [k/2 for k in h]
#       
# # error due to h
# size = [3100,2200,1300,1040]
# e = []
# for dt,s in zip( h,size):
#     E = diffusion(dt,s)
#     e.append(E)
# # error du to h/2
# e_hf = []
# for dt,s in zip( h_hf,size):
#     E = diffusion(dt,s)
#     e_hf.append(dt)
#           
# P = []
# P_hf = []
# Eh = []
#     
# # convergence rate calculation
# for i in range(len(e)):
#     p_hf = abs(log(e_hf[i]/e_hf[i-1])/log(h_hf[i-1]/h_hf[i]))
#     P.append(p)
#     P_hf.append(p_hf)
#      
#  
dt = 0.045/2
s = 3100
# 
diffusion(dt,s)
#print e_hf
# 
#print"convergence rate using h/2:",P_hf
#print "Error", e_hf
#   
  


    
    





