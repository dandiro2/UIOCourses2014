
from matplotlib import rc
#rc('text', usetex=True)

from math import pi, sin
import numpy as np
import sympy as sp
import matplotlib.pyplot as pl
from scipy.integrate import quad


   
# define the basis functions phi0,phi1,and phi2 
def phi0(x):
    if 0<= x <= pi/2.:
        return ( (pi/2.) -x)/(pi/2.)
    else:
        return 0
    
    
def phi1(x):
    if 0<= x <= pi/2.:
        return (x)/(pi/2.)
    elif pi/2 <= x <=pi:
        return (pi-x)/(pi/2.)
    
    else:
        return 0
    
def phi2(x):
    if pi/2.<= x <= pi:
        return (x-(pi/2.))/(pi/2.)
    else:
        return 0
    
def sin_approx_P1():
    """
    Aprroximate sin(x) with P1 element od size pi/2
    on [0,pi]
    """ 
    x = sp.Symbol("x")
    A = np.zeros((3,3))
    b = np.zeros(3)
    
    
    r = pi/2
    A[0,0] = sp.integrate(   ((r-x)/(r))*((r-x)/(r)) ,(x,0,r)  )
    A[1,1] = sp.integrate(   ((x)/(r))*((x)/(r)) ,(x,0,r)  ) + sp.integrate(   ((pi-x)/(r))*((pi-x)/(r)) ,(x,r,pi)  )
    A[2,2] = sp.integrate(   ((-r+x)/(r))*((-r+x)/(r)) ,(x,r,pi)  )
    
    A[2,1] = sp.integrate(   ((pi-x)/(r))*((x-r)/(r)) ,(x,r,pi)  )
    
    A[0,1] = sp.integrate(   ((r-x)/(r))*((x)/(r)) ,(x,0,r)  )
    
    A[1,2] = A[2,1]
    A[1,0] = A[0,1]
    
    
    
    b[0] = sp.integrate(   ((r-x)/(r))*(sp.sin(x)) ,(x,0,r)  )
    
    b[1] = sp.integrate(   ((x)/(r))*(sp.sin(x)) ,(x,0,r)  ) + sp.integrate(   ((pi-x)/(r))*(sp.sin(x)) ,(x,r,pi)  )
    
    b[2] = sp.integrate(   ((x-r)/(r))*(sp.sin(x)) ,(x,r,pi)  )

    c = np.linalg.solve(A,b)
    
    z = np.linspace(0,pi,1000)
    h = c[0]*(   (r-z)/r )  +c[1]*(z/r)  + c[2]*((z-r)/r)
    h1 = c[0]*(   (r-z)/r )  +c[1]*((pi-z)/r)  + c[2]*((z-r)/r)

    print c
    print ""
    print A
    
    print b


sin_approx_P1()
    
    
    
# omega0 = np.linspace(0,pi/2,1000)
# omega1 = np.linspace(pi/2.,pi,1000)
# Phi0 = []
# Phi11 = []
# Phi12 = []
# Phi2 = []
# for x in omega0:
#     Phi0.append(phi0(x))
#     
# for x in omega0:
#     Phi11.append(phi1(x))
# for x in omega1:
#     Phi12.append(phi1(x))
# for x in omega1:
#     Phi2.append(phi2(x))
# pl.plot(omega0, Phi0, omega0,Phi11,"r",omega1,Phi12,"r", omega1,Phi2)
# pl.title("basis functions $\phi_{0}$,  $\phi_{1}$, $\phi_{2}$")
# pl.show()

