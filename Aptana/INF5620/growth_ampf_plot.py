
from sympy import *
import numpy as np
from matplotlib import *
import pylab as pl
from math import exp



def amplificationFactor():
    
    theta, a,dt,p = symbols("theta a dt p")
    A = ( 1-(1-theta)*a*dt)/( 1+theta*a*dt )
    
    ampf = A
    half = Rational(1,2)
    ampff = A.subs({a*dt:p,theta:0})
    print ampff
    ampff = lambdify(p,ampff)
    
    ampfb = A.subs({a*dt:p,theta:1})
    ampfb = lambdify(p,ampfb)
    
    ampfcn = A.subs({a*dt:p,theta:0.5})
    ampfcn = lambdify(p,ampfcn)
    
    
#     ampff  = lambdify(dt,ampf.subs({theta:0,a*dt:p}))
#     ampfb  = lambdify(dt,ampf.subs({theta:1,a*dt:p})   )
#     ampfcn = lambdify(dt,ampf.subs({theta:half,a*dt:p} )  )

    
    
    p = np.linspace(-3,-1,100)
    pb1 = np.linspace(-1.5,-0.99,100)
    pb2 = np.linspace(-1.1,-0.5,100)
    
    pcn1 = np.linspace(-3,-1.99,100)
    pcn2 = np.linspace(-2.1,-0.5,100)

    y = 1-p
    
    Ae = np.exp(-p)
    pl.subplot(221)
    pl.plot(p,Ae)
    pl.title("exact")
    pl.grid(True)
    
    pl.subplot(222)
    pl.plot(p,ampff(p))
    pl.grid(True)
    pl.title("forward")
    
    pl.subplot(223)
    pl.plot(pb1,ampfb(pb1),"r--",pb2,ampfb(pb2), "r--" )
    pl.grid(True)
    pl.title("backward")
    
    pl.subplot(224)
    pl.plot(pcn1,ampfcn(pcn1),"r--",pcn2,ampfcn(pcn2), "r--" )
    pl.grid(True)
    pl.title("backward")
    
    pl.show()
    
    
#     pl.plot(dt,Ae,"b-", dtf,ampff(dtf),"r-",dtb1,ampfb(dtb1),"o",dtb2,ampfb(dtb2),"o",dtcn1,ampfcn(dtcn1),"r--o",dtcn2,ampfcn(dtcn2),"r--o"  )
#     pl.legend(["ecact","forward","backward","crank nicolson"])
#     pl.grid(True)
#     
#     pl.show()
    
    
    
    
#     print " ampf forrward: ", ampff
#     print " ampf backward: ", ampfb
#     print " ampf cN :" ,ampfcn

if __name__=="__main__":
    print "ok"
    
    
