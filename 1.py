from dolfin import *
import numpy as np

h = 8
tau = 1.0/128.0

mesh = UnitSquareMesh(h,h)
V = FunctionSpace(mesh,"CG",1)
u = TrialFunction(V)
v = TestFunction(V)


k = Constant(1.0)
f = Constant(1.0)
T = 1.0

a = (1/tau)*u*v*dx + inner(k*grad(u),grad(v))*dx

u0 = Function(V)
f = Constant(1)

#define linear part
L = (u0 / tau)*v*dx+ f*v*dx

u = Function(V)

#boundary
def Bound(x,on_boundary):
	return on_boundary

DBC = DirichletBC(V,Constant(0),Bound)

#file = File("1.pvd")
ar_max = np.empty(128)
i = 0
for t in np.arange(0.0,T,tau):
	solve(a == L, u, DBC)
	assign(u0, u)
	ar_max[i] = u0.vector().max()
	i=i+1
	#file << (u, t)

print ar_max

