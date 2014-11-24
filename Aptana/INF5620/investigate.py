from wavesolver import *
from scipy.integrate import dblquad
import numpy as np
# Initial Gauss, (with y-dependency)
def I(x,y,I_0=0.8,I_a=0.4,I_m=0.9,I_s=0.2): 
    return I_0 + I_a*np.exp(-((x-I_m)/I_s)**2-((y-I_m/2)/I_s)**2)

# Bottom_Gauss
def B_G(x,y,B_0=0.0,B_a=0.6,B_mx=0.2,B_my=0.5,B_s=0.1,b=4):
    return B_0 + B_a*np.exp(-((x-B_mx)/B_s)**2 - ((y-B_my)/(b*B_s))**2)

# Bottom cos
def B_c(x,y,B_0=0.6,B_a=0.1,B_mx=0.1,B_my=0.5,B_s=0.6):
    A = np.sqrt(x**2+y**2) <= B_s
    return A*(B_0 + B_a*np.cos(np.pi*(x-B_mx)/(2*B_s))*np.cos(np.pi*(y-B_my)/(2*B_s)))

# bottom hill
def hill(x,y):
    return (y-0.2)*(y>0.2)    # hill function
    
bottom = lambda x,y:0    # set sea-bottom profile

Lx = 1.0
Ly = 1.0
H_eq,err = dblquad(I,0,Ly,lambda x: 0,lambda x:Lx) # compute still water heigth
H_eq=H_eq/(Lx*Ly)

I = lambda x,y:1
g = 9.81
def f(x,y,t):
    Cg = np.sqrt(g*(H_eq-bottom(x,y)))
    A = 80
    yline = (y>0.49)*(y<0.51)
    yr = (y>0.5)*(y<0.52)
    yl = (y>0.48)*(y<0.5)
    pos = abs(x-Cg*t)<0.1
    pos2 = abs(x+0.1-Cg*t)<0.1
    return A*0.2*pos*yline-A*pos2*0.2*yline+2*A*np.sin(300*t)*yr*pos2+2*A*np.sin(300*t)*yl*pos2

#0*y-60*t)*((x-0.5)/(t+0.001)>0.1)
b = 0
q = lambda x,y: g*(H_eq-bottom(x,y))
V = lambda x,y:0
#wave = wavesolver(I,f,q,V,b)
solve(I,f,q,b)