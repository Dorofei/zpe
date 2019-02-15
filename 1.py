from dolfin import *
import numpy as np


mesh = UnitSquareMesh(10,10)
V = FunctionSpace(mesh,"CG",1)
u = TrialFunction(V)
v = TestFunction(V)

tau = 0.01;
k = Constant(1.0)
f = Constant(1.0)
T = 1.0

a = 1/tau*u*v*dx + inner(k*grad(u),grad(v))*dx

u0 = Function(V)
f = Constant(1)

#define linear part
L = (u0 / tau)*v*dx+ f*v*dx

u = Function(V)

#boundary
def Bound(x,on_boundary):
	return on_boundary

DBC = DirichletBC(V,Constant(0),Bound)

file = File("./result1/1.pvd")
for t in np.arange(0.0,T,tau):
	solve(a == L, u, DBC)
	assign(u0, u)
	file << (u, t)

