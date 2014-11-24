from dolfin import *

#define mesh
mesh = UnitIntervalMesh(2)
#mesh = IntervalMesh(20,0,1)
V = FunctionSpace(mesh,'Lagrange',2)

# set Dirichlet boundary condition
u0 = Expression('1')
def u0_boundary(x,on_boundary):
    return on_boundary

bc = DirichletBC(V,u0,u0_boundary)

# varional formulation
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(1)
a = -inner(grad(u),grad(v))*dx
L = f*v*dx

# solution of the variatiobal problem
u = Function(V)
solve(a==L,u,bc)

#plot solution
plot(u)

# view solution
interactive()
