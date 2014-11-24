from dolfin import *
import numpy as np
from math import log
import pylab as pl
def diffusion(I,f,dt,T,t,rho,size,degree,d,dimension,alpha = None):
    
    # Create mesh and define function space
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    elif dimension == 2:
        mesh == UnitSquareMesh(size,size)
    elif dimension == 3:
        mesh = UnitCube(size,size,size)
    
    # function space
    V = FunctionSpace(mesh, "Lagrange", d)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_1 = interpolate(I,V) 
    if alpha == None:
        
        def alpha(u):
            return 1
    
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

    return  u
    

 
  
def ConstantSolutionVerification(d):
    """
    output:
    ##########################################################
    The L2 norm for the constant solution:  9.55040411904e-15
    ##########################################################

    """
    degree, size, dimension = 1,2000,1
    mesh = UnitIntervalMesh(size)
    V = FunctionSpace(mesh,'Lagrange',d)
    C = 4.5
    I = Constant(C)
    f = Constant(0.0)
    rho, dt, T, t = 1000, 0.01, 0.1, 0
    uexact = Constant(C)
    ua = diffusion(I,f,dt,T,t,rho,size,degree,d,dimension)
    
    U = interpolate(ua,V)
    W = interpolate(uexact,V)
    e = np.abs(W.vector().array()-U.vector().array()).max()
      
    L_2norm = np.sqrt(np.sum(e**2)/ua.vector().array().size)
    tol = 1e-14
    
    assert L_2norm <= tol
    print""
    print"##########################################################"
    print"The L2 norm for the constant solution: ", L_2norm
    print"##########################################################"


def nonlinearVerification(d):
    """
    ##################################################################
    The L2 norm for the nonlinear solution :  2.85931976232e-15
    ###################################################################

    """
    rho, dt, T, t = 1, 0.01, 0.1, 0
    def alpha(u):
        return 1+u*u
    degree, size, dimension = 1,2000,1
    mesh = UnitIntervalMesh(size)
    V = FunctionSpace(mesh,'Lagrange',1)
    I = Expression('t*x[0]*x[0]*(0.5-(1./3)*x[0])',t=0)
    f = Expression('rho*pow(x[0],2)*(-2*x[0] + 3)/6 - \
        (-12*t*x[0] + 3*t*(-2*x[0] + 3))*(pow(x[0],4)*pow((-dt + t),2)*pow((-2*x[0] + 3),2) + 36)/324 - \
        (-6*t*pow(x[0],2) + 6*t*x[0]*(-2*x[0] + 3))*(36*pow(x[0],4)*pow((-dt + t),2)*(2*x[0] - 3) + \
        36*pow(x[0],3)*pow((-dt + t),2)*pow((-2*x[0] + 3),2))/5832',
        rho=rho, t=t, dt=dt)
    
    uexact = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)', t=0.1)
    ua = diffusion(I,f,dt,T,t,rho,size,degree,d,dimension,alpha)
    U = interpolate(ua,V)
    W = interpolate(uexact,V)
    e = np.abs(W.vector().array()-U.vector().array()).max()
      
    L_2norm = np.sqrt(np.sum(e**2)/ua.vector().array().size)
    tol = 1e-14
    
    assert L_2norm <= tol
    print""
    print"##################################################################"
    print"The L2 norm for the nonlinear solution : ", L_2norm
    print"###################################################################"

def analitycalVerification():
    """
    output:
    ###################################################################################################
    ratio E/h:  [1.4823146933763534, 1.483393109737996, 1.4789673956393048, 1.4805356573640203]
    convergence rate : [0.9999999999999999, 0.9999999999999994, 0.999999999999999, 1.0]
    #######################################################################################################

    """
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
    print"###################################################################################################"  
    print"ratio E/h: ",Eh 
    print"convergence rate :",P_hf
    print"#######################################################################################################" 
      
#analitycalVerification()
nonlinearVerification(1) 
#ConstantSolutionVerification(1)




