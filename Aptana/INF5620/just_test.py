# NB: the function must be called giving degree d and mesh sizes
# as parameters in the terminal,
# e.g. python fem_project.py 1 100 100
# for a 1D problem with the unit interval divided in 100 intervals

from dolfin import*
import sys 
import numpy as np
import pylab as pl

def diffusion_solver(u_exact, I, f, alpha, rho=1, dt=0.01, T=0.1, divisions=None):
    # Define from the terminal the polynomial degree and
    # the type of domain (1D, 2D or 3D), saved in number of divisions. 
    """
    degree = int(sys.argv[1])
    if divisions == None:
        divisions = [int(arg) for arg in sys.argv[2:]]
    d = len(divisions)
    domain_type = [UnitIntervalMesh, UnitSquareMesh, UnitCubeMesh]    
     """
    #create mesh
    #mesh = domain_type[d-1](*divisions)
    mesh = UnitIntervalMesh(12)
    
    #create function space
    V = FunctionSpace(mesh, 'Lagrange', 1)
    
    #load the first "previous known" u, which will be the approximation of 
    #the last time step in order to start the Picard iterations
    u_1 = interpolate(I,V)
    
    #define the timestepping
    t = dt

    while t <= T:
        #evaluate the function at time t
        f.t = t
        
        #define trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)
        
        #the variational formulation for Picard method
        a = u*v*dx +dt/rho*inner(alpha(u_1)*grad(u),grad(v))*dx
        L = (u_1 + f*dt/rho)*v*dx

        u = Function(V)
        solve(a == L, u)
        
        #update values for next iteration
        t += dt
        u_1.assign(u)

    #visualize
    #plot(u)
    u_e = interpolate(u_exact,V)
    #plot(u_e)

#     U = u.vector().array()
#     w = u_e.vector().array()
#     x = np.linspace(0,1,len(U))
#     pl.plot(x,U,'-',x,w,'ro')
#     pl.show()
    
    #interactive()

    #verify the approximation with the maximum error
    #e_max = np.abs(u_e.vector().array() - u.vector().array()).max()
    #print 'Max error, t=%.2f: %-10.3f' % (t, e_max)
       
    #compute the discrete L2 norm of the solution at the nodes, for the last time step
    e = u_e.vector().array() - u.vector().array()      
    L_2norm = np.sqrt(np.sum(e**2)/u.vector().array().size)
    return L_2norm, t


if __name__ == '__main__':    
    
    def constant():
        C = 5.0
        I = Constant(C)
        f = Constant(0.0)
        alpha = Constant(1.0)

        u_exact = Constant(C)
        diffusion_solver(u_exact, I, f, alpha)
        
    
    def analytical():
        I = Expression('cos(pi*x[0])')
        f = Constant(0.0)
        alpha = lambda u: 1.0

        u_exact = Expression('exp(-pow(pi,2)*t)*cos(pi*x[0])', t=0)
         
        h = []  
        E = []
        t = []
        for N in [4, 8, 16, 32]:
            divisions = [N, N]
            dt = (1.0/N)**2
            print '###########Solving for N = %d with dt = %f###########' % (N, dt)
            h.append(dt)
            error, time = diffusion_solver(u_exact, I, f, alpha, dt=dt, divisions=divisions)
            E.append(error)
            t.append(time)

        print '#########################################################'
        for i in range(len(E)):
            n = E[i]/h[i]
            print 'At time %f: E = %f h = %f E/h = %f' % (t[i], E[i], h[i], n)


    def manufactured():
        I = Constant(0.0)
        rho = 1.0
        f = Expression('-rho*pow(x[0],3)/3 + rho*pow(x[0],2)/2 \
        + 8*pow(t,3)*pow(x[0],7)/9 - 28*pow(t,3)*pow(x[0],6)/9 \
        + 7*pow(t,3)*pow(x[0],5)/2 - 5*pow(t,3)*pow(x[0],4)/4 \
        + 2*t*x[0] - t', rho=rho, t=0)
        def alpha(u):
            return 1 + u**2
        #alpha = lambda u: 1 + u**2
        u_exact = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)', t=0.1)
        diffusion_solver(u_exact, I, f, alpha)


    def verification_of_manufactured():        
        I = Constant(0.0)
        rho = 1.0
        f = Expression('rho*pow(x[0],2)*(-2*x[0] + 3)/6 - \
        (-12*t*x[0] + 3*t*(-2*x[0] + 3))*(pow(x[0],4)*pow((-dt + t),2)\
        *pow((-2*x[0] + 3),2) + 36)/324 \
        - (-6*t*pow(x[0],2) + 6*t*x[0]*(-2*x[0] + 3))*(36*pow(x[0],4)\
        *pow((-dt + t),2)*(2*x[0] - 3) \
        + 36*pow(x[0],3)*pow((-dt + t),2)*pow((-2*x[0] + 3),2))/5832', 
        rho=rho, t=0, dt=0.01)
        alpha = lambda u: 1 + u**2
        u_exact = Expression('t*pow(x[0],2)*(0.5 - x[0]/3.)', t=0.1)

        error, time = diffusion_solver(u_exact, I, f, alpha)
        print '#########################################################'
        dt = 0.01
        n = error/dt
        #print 'At time %f: E = %f h = %f E/h = %f' % (time, error, dt, n)
        


    def gaussian():
        rho = 1.0
        dt = 0.01
        T = 0.1
        sigma = 1.0
        beta = 0.3

        I = Expression('exp(-1/(2*pow(sigma,2))*(pow(x[0],2) + pow(x[1],2)))', sigma=sigma)
        f = Constant(0.0)
        alpha = lambda u: 1 + beta*u**2
        u_exact = I

        diffusion_solver(u_exact, I, f, alpha)

    manufactured()
    verification_of_manufactured()
