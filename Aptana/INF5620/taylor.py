from sympy import *
import numpy as np
from sympy.mpmath import *
import pylab as pl

def f(x):
    return exp(-x)

q = 100
z = np.linspace(0,4,q)
x0  = 2
deg = 6
p = taylor(f,x0,deg)

App = []
for x in z:
    r = polyval(p[::-1], x-x0)
    App.append(r)

pl.plot(z,App,"bo", z,np.exp(-z))
pl.show()

