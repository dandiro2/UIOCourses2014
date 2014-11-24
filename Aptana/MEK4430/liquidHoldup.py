"""
solve the holdup equation for sigma

"""
import numpy as np
import pylab as pl
from math import pi, sin, cos
from scipy.optimize import fsolve, bisect

def liquidHoldup(theta,Usg,Usl,beta,string=None):
    """
    this function takes as input the parameters theta 
    the superficial velocity of liquid and gas: Usl, Usg 
    and the string "print".
    
    usage: 
    to print results, enter the string "print":
    
    thata = 1; Usg = 3; Usl = 0.3; beta = 0
    liquidHoldup(thata,Usg,Usl,beta,"print")
    
    output:
    ################################################
    parameters entered:
    superficial liquid velocity Usl =  0.3
    superficial gas velocity Usg    =  3
    angle beta                      =  0
    #################################################
    computation results:
    liquid holdup             = 0.837510724346
    gas holdup                = 0.162489275654
    
    in-situ liquid velocity   = 0.358204368349 m/s
    in-situ gas velocity      = 18.4627569292 m/s
    ################################################
    
    to return the value of sigma to be used later in a computation:
    do not specify the string print:
    
    usage:
    sigma = liquidHoldup(thata,Usg,Usl,beta)

    """
    
    def holdupEquation(sigma):
        """
        This function return the holdup equation as a function of
        sigma. The holdup equation is used as an input
        by the fsolve() python function to return sigma.
        the python function fsolve() solves nonlinear equation.
        """
        #physical parameter values
        rhoL, rhoG = 1000, 50       # density

        D          = 0.1              # pipe diameter
        muL, muG   = 0.001, 0.00001 # dynamic viscosity
    
        # in-situ velocity
        alphaL = (sigma-sin(sigma))/(2*pi)
        UL = Usl/(alphaL)
        UG = Usg/(1-alphaL)
        
        # physical boundary of the pipe
        A  = (pi*D**2)/4                    # cross sectional area of pip
        AL = (A/2*pi)*(sigma-sin(sigma))    # cross sectional area of liquid
        AG = A-AL                           # cross sectional area of gas
        
        SL = (D*sigma)/2  
        SG = D*(pi-0.5*sigma)
        SI = D*sin(0.5*sigma)
        
        DL = (4*AL)/SL
        DG = (4*AG)/(SI+SG)               
        
        #shear stress
        ReG = (rhoG*UG*DG)/muG # gas Reynold number
        ReL = (rhoL*UL*DL)/muL # liquid Reynold number
    
        #blasius friction factor
        fG = 0.316/(ReG**0.25)            # gas friction factor
        fL = 0.316/(ReL**0.25)            # liquid friction factor
    
        #shear stress
        tauG = (1./8)*rhoG*fG*UG**2         # gas shear sress
        tauL = (1./8)*rhoL*fL*UL**2         # liquid shear stress
        tauI = (1./8)*rhoG*theta*fG*(UL-UG)**2    # interfacial shear stress
        g = 9.81 # gavitational constant
        
        #holdup equation as a function of sigma
        equation = ( ((tauG*SG)/AG)-((tauL*SL)/AL)   ) + tauI*SI*((1./AG)+(1./AL) )+g*sin(beta)*((rhoG/AG)-(rhoL/AL) )
        return equation
    

    sig = bisect(holdupEquation,0.1,pi/2.)
    #sig = fsolve(holdupEquation,0.8)
     
    #compute liquid holdup
    alphaL = (sig-sin(sig))/2.*pi
     
    #compute in-situ velocities
    UL = Usl/alphaL
    UG = Usg/(1-alphaL)
     
    #print result if string print is an iput to function liquidHolup()
    if string !=None:
        print"################################################"
        print"parameters entered:"
        print"superficial liquid velocity Usl = ", Usl
        print"superficial gas velocity Usg    = ", Usg
        print"angle beta                      = ", beta
        print"theta                           = ", theta
        print"#################################################"
        print"computation results:"
        print"liquid holdup             =", alphaL
        print"gas holdup                =", 1-alphaL
        print""
        print"in-situ liquid velocity   =",UL, "m/s"
        print"in-situ gas velocity      =",UG, "m/s"
        print"################################################"
        print"h", (0.1/2)*(1-cos(sig/2))
         
         
    
    return sig
    

def pressureGradient():
    """
    given the sigma compute the pressure gradient
    """
    theta, Usg, Usl, beta = 1,3,0.3,0
    sigma = liquidHoldup(theta,Usg,Usl,beta)
    rhoL, rhoG = 1000, 50       # density

    D          = 0.1              # pipe diameter
    muL, muG   = 0.001, 0.00001 # dynamic viscosity
    
    # in-situ velocity
    alphaL = (sigma-sin(sigma))/(2*pi)
    UL = Usl/(alphaL)
    UG = Usg/(1-alphaL)
        
    # physical boundary of the pipe
    A  = (pi*D**2)/4                    # cross sectional area of pip
    AL = (A/2*pi)*(sigma-sin(sigma))    # cross sectional area of liquid
    AG = A-AL                           # cross sectional area of gas
        
    SL = (D*sigma)/2  
    SG = D*(pi-0.5*sigma)
    SI = D*sin(0.5*sigma)
        
    DL = (4*AL)/SL
    DG = (4*AG)/(SI+SG)               
        
    #shear stress
    ReG = rhoG*UG*DG/muG # gas Reynold number
    ReL = rhoL*UL*DL/muL # liquid Reynold number
    print ReG,ReL
    
    fG = 0.316/(ReG**0.25)            # gas friction factor
    fL = 0.316/(ReL**0.25)            # liquid friction factor
    
    tauG = (1./8)*rhoG*fG*UG**2         # gas shear sress
    tauL = (1./8)*rhoL*fL*UL**2         # liquid shear stress
    tauI = (1./8)*rhoG*theta*fG*(UL-UG)**2    # interfacial shear stress
    g = 9.81 # gavitational constant
    
    PLI = (1./AL)*(tauL*SL-tauI*SI) +(rhoL*g*sin(beta))/(AL)
        
    PGI = (1./AG)*(tauG*SG+tauI*SI) +(rhoG*g*sin(beta))/(AG)
    
    print"parameters used for computation:"
    print"theta = ", theta
    print"Usl   = ", Usl
    print"Usg   = ", Usg
    print"beta  = ", beta
    print"##################################################"
    print"pressure gradient for liquid dPLi/dx", PLI
    print"pressure gradient for Gas dPGi/dx   ", PGI
    print"###################################################"
    print"DL",DL
    print"DG", DG

        
     
    
#pressureGradient()
theta, Usg,Usl,beta = 1,3,0.3,0
liquidHoldup(theta,Usg,Usl,beta,"print")




    
    
    
    