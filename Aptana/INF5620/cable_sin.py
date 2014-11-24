import numpy as np
import sympy as sp
from scipy import integrate
import pylab as pl

def cable_sin(N,str):
    """
    solve the differential equation :
    u'' = 1;  with two different boundary conditions:
    part 1:
    u(0) = 0; u'(1) = 0 (Dirichlet bondary condition at 0 and Newmann boundary condition at 1.
    basis functions: phi_i(x) = sin( (i+1)*pi*x*0.5)
    u_exact(x) = 0.5*x**2-x
    
    part 2)
    2)u(0) = u(1) = 1  (Dirichlet boundary condition at 0 and 1
    basis functions: phi_i(x) = sin( (i+1)*pi*x)
    analytical solution :
    u(x) = 0.5*x**2-0.5x
    
    Usage of the function:
    
    cable_sin(0,"part1")   for part 1 out put is :
    
    ######################################################################
    error in maximum deflection at x = 1 for N = 0 :  0.0160245509312

    exact solution at x = 1:      -0.5

    numerical solution at x = 1:  -0.516024550931
    #######################################################################
    
    
    cable_sin(0,"part2")   for part 2 out put is :
    #########################################################
    solution of symetry problem at u(1):       -0.516024550931

    solution of non-symetry problem at u(1/2): -0.516024550931
    #########################################################
    """
    x = sp.Symbol("x")

    
    
    # basis function and their derivatives for part 1:
    if str == "part 1":
        phi = [sp.sin((i+1)*sp.pi*0.5*x) for i in range(N+1)]
        phidiff = [sp.diff(phi[i],x) for i in range(N+1)]
        
        #stress matrix A and load matrix b
        A = np.zeros((N+1,N+1))
        b = np.zeros(N+1)
        
        for i in range(N+1):
            integrand1 = sp.lambdify(x,phi[i])
            b[i],error = integrate.quad(integrand1 ,0,1)  
             
            for j in range(N+1):
                integrand2 = sp.lambdify(x,phidiff[i]*phidiff[j])
                A[i,j],error = integrate.quad(integrand2 , 0,1)     
        
        # coefficient matrix c    
        c = np.linalg.solve(-A,b)
    
    
        # finite elemtn solution u
        U = sum( [c[i]*phi[i] for i in range(N+1) ])
        u = sp.lambdify(x,U,"numpy")
    
        # analytical solution
    
    
        def ue(x):
            return 0.5*x**2 - x
        print "######################################################################"
        print "N = ", N
        print "error in maximum deflection at x = 1 : ", abs( ue(1)-u(1) )
        print ""
        print "exact solution at x = 1:     ",ue(1)
        print ""
        print "numerical solution at x = 1: ",u(1)
        print "#######################################################################"
        
        C = [ abs(c[j]/c[j-1]) for j in range(1,N+1)]
        print C
        Z = np.linspace(0,15,N+1)
        pl.plot(C)
        pl.xlabel("N")
        pl.ylabel("cj/cj-1")
        pl.title("variation of cj/cj-1 with N")
        pl.show()
        
        
           
    
    elif str == "part 2":
        # basis functions and there derivatives for part 2
        
        phi1 = [sp.sin((i+1)*sp.pi*0.5*x) for i in range(N+1)]
        phidiff1 = [sp.diff(phi1[i],x) for i in range(N+1)]
        
        phi2 = [sp.sin((i+1)*sp.pi*x) for i in range(N+1)]
        phidiff2 = [sp.diff(phi2[i],x) for i in range(N+1)]
        
        A1 = np.zeros((N+1,N+1))
        b1 = np.zeros(N+1)
        A2 = np.zeros((N+1,N+1))
        b2 = np.zeros(N+1)
        
        for i in range(N+1):
            integrand1 = sp.lambdify(x,phi1[i])
            integrand2 = sp.lambdify(x,phi2[i])

            b1[i],error = integrate.quad(integrand1 ,0,1)  
            b2[i],error = integrate.quad(integrand2 ,0,1)  
   
        for j in range(N+1):
            integrand11 = sp.lambdify(x,phidiff1[i]*phidiff1[j])
            A1[i,j],error = integrate.quad(integrand11 , 0,1)     
            
            integrand22 = sp.lambdify(x,phidiff2[i]*phidiff2[j])
            A2[i,j],error = integrate.quad(integrand22 , 0,1)     
    
        c1 = np.linalg.solve(-A1,b1)
        c2 = np.linalg.solve(-A2,b2)
        print c1[0]
        print c2[0]
        print c1[0]/c2[0]
        from math import sqrt
        print 2./sqrt(2)
        U1 = sum( [c1[i]*phi1[i] for i in range(N+1) ])
        u1 = sp.lambdify(x,U1,"numpy")
        
        U2 = sum( [c1[i]*phi2[i] for i in range(N+1) ])
        u2 = sp.lambdify(x,U2,"numpy")
        
        print"#########################################################"
        print"solution of symetry problem at u(1):      "  ,u1(1)
        print""                                                     
        print"solution of non-symetry problem at u(1/2):" ,u2(0.5)
        print"#########################################################"
    

        
        
    

    

        
        
cable_sin(0,"part 2")




