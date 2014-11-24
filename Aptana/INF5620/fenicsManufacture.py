from dolfin import *
import numpy as np
import pylab as pl

def diffusion(dimension):
    
    # Create mesh and define function space
    size = 1200
    if dimension == 1:
        mesh = UnitIntervalMesh(size)
    elif dimension == 2:
        mesh == UnitSquareMesh(size,size)
    elif dimension == 3:
        mesh = UnitCube(size,size,size)
    
    # function space
    I = Constant(0.0)
    rho = 1.0
    d = 1
    f = Expression('-rho*pow(x[0],3)/3 + rho*pow(x[0],2)/2 \
    + 8*pow(t,3)*pow(x[0],7)/9 - 28*pow(t,3)*pow(x[0],6)/9 \
    + 7*pow(t,3)*pow(x[0],5)/2 - 5*pow(t,3)*pow(x[0],4)/4 \
    + 2*t*x[0] - t', rho=rho, t=0)
    def alpha(u):
        return 1 + u**2
        #alpha = lambda u: 1 + u**2
    u_exact = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)', t=2)
    V = FunctionSpace(mesh, "Lagrange", d)
    u = TrialFunction(V)
    v = TestFunction(V)
    u_1 = interpolate(I,V) 
    
    T = 2
    dt = 0.01
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

    ue = interpolate(u_exact,V)
#     U = u.vector().array()
#     W = ue.vector().array()
#     x = np.linspace(0,1,len(U))
#     
#     pl.plot(x,U,'-',x,W,'ro')
#     pl.legend(['fenics','manufacture'])
#     pl.show()
    plot(ue)
    plot(u)
    interactive()

dimension = 1
diffusion(dimension)
