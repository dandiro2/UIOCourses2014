from dolfin import *
import numpy as np
from math import log

def analitycalVerfication():
    def diffusion(dt,size):    
        # Create mesh and define function space
        mesh = UnitIntervalMesh(size)
        V = FunctionSpace(mesh, "Lagrange", 1)
        R = FunctionSpace(mesh, "R", 0)
        W = V * R
        # Define variational problem
        (u, c) = TrialFunction(W)
        (v, d) = TestFunctions(W)
        f = Constant(0.0)
        g = Constant(0.0)
        T = 0.05
        rho = 1.
        I = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0)
        u_1 = interpolate(I,V) 
        alpha = Constant(1.0)
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
        
        #plot(u)
        
        u_e = interpolate(I, V)
        #plot(u_e)
        #plot(u,interactive=True)
        #plot(u_e,interactive=True)
        e = np.abs(u_e.vector().array()-u.vector().array()).max()
          
        L_2norm = np.sqrt(np.sum(e**2)/u.vector().array().size)
        return e
    
    
    #
    h = [0.00045,0.00040,0.00035,0.00030]
    h_hf = [k/2 for k in h]
      
    # error due to h
    size = [16,11,7,6]
    e = []
    for dt,s in zip( h,size):
        E = diffusion(dt,s)
        e.append(E)
    # error du to h/2
    e_hf = []
    for dt,s in zip( h_hf,size):
        E = diffusion(dt,s)
        e_hf.append(dt)
          
    P = []
    P_hf = []
    Eh = []
    
    # convergence rate calculation
    for i in range(len(e)):
        p_hf = abs(log(e_hf[i]/e_hf[i-1])/log(h_hf[i-1]/h_hf[i]))
        P.append(p)
        P_hf.append(p_hf)
      
    # calculation of ratio E/h : need to vary the mesh zise as well
    
    ee = []
    size = [3100,2200,1300,1040]
    for dt,s in zip(h,size):
        E = diffusion(dt,s)
        ee.append(E)
    
    for i,j in zip(ee,h):
        eh = i/j
        Eh.append(eh)
    print""
    print"#########################################"  
    print " ratio E/h: ",Eh 
    print"convergence rate :",P_hf
    print"########################################" 
      


    
    





