"""
solve for Cj by applying Galerkin mathod 
to the problem:
u'' = 1 if x < 0.5
u'' = 0 if x > 0.1
with boundary condition:
u(0) = u(1) = 0

output:
  raise LinAlgError("Singular matrix")
numpy.linalg.linalg.LinAlgError: Singular matrix
"""



import sympy as sp
import numpy as np
from scipy.integrate import quad
N = 10
A = np.zeros((N+1,N+1))
A1= np.zeros((N+1,N+1))
b = np.zeros(N+1)
b1= np.zeros(N+1)
x = sp.Symbol('x')
phi      = [ sp.sin((i+1)*sp.pi*x) for i in range(N+1)]
phi_diff2 = [sp.diff(sp.diff(phi[i],x),x) for i in range(N+1)]

for i in range(N+1):
    def integrand1(x):
        return 1
    b[i], error = quad(integrand1,0,0.5)
    for j in range(N+1):
        integrand2 = sp.lambdify(x,phi_diff2[i])
        A[i,j],error = quad(integrand2 , 0,1)
    
print A
c1 = np.linalg.solve(A,b)

        
