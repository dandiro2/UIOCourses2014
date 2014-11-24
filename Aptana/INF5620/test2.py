from dolfin import *

# Create mesh and define function space
#mesh = UnitSquareMesh(64, 64)
mesh = UnitIntervalMesh(60)
V = FunctionSpace(mesh, "Lagrange", 1)
R = FunctionSpace(mesh, "R", 0)
W = V * R

# Define variational problem
(u, c) = TrialFunction(W)
(v, d) = TestFunctions(W)
#f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)")
#g = Expression("-sin(5*x[0])")
T = 1
dt = 0.05
rho = 1.
C = 1.6
f = Constant(0.0)
g = Constant(0.0)
I = Constant(C)
alpha = interpolate(I,V)
u_1 = interpolate(I,V)
a = ((dt/rho)*inner(alpha*grad(u), grad(v)) +u*v+ (dt/rho)*c*v + u*d)*dx
L = (dt/rho)*f*v*dx + g*v*ds +u_1*v*dx
A = assemble(a)
b = assemble(L)
# Compute solution
w = Function(W)
#u_1 = Function(V)
t= dt
while t<=T:
    solve(A,w.vector(),b)
    (u, c) = w.split(True)
    u_1.assign(u)
    t+=dt

u_ex = Expression('exp(-pi*pi*t)*cos(pi*x[0])',t=0.05)
ue = interpolate(u_ex,V)
# Plot solution
plot(u, interactive=True)
