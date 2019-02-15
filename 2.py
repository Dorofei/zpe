from dolfin import *
import numpy as np
set_log_level(WARNING)

t_array = np.array([64,128,256,512]) 	#1/64, 1/128, 1/256, 1/512
h_array = np.array([8,16,32,64])		#1/8, 1/16, 1/32, 1/64
T = 1.0/8.0

for h in np.nditer(h_array):
	for tt in np.nditer(t_array):
		print "h = 1 /",h,"tt = 1 / ",tt

		tau = 1.0/tt
		mesh = UnitSquareMesh(int(h),int(h))
		V = FunctionSpace(mesh,"CG",1)
		u = TrialFunction(V)
		v = TestFunction(V)


		k = Constant(1.0)
		f = Constant(1.0)


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

		file = File("./result2/"+str(h)+"-"+str(tt)+".pvd")
		ar_max = np.empty(tt)
		i = 0
		for t in np.arange(0.0,T,tau):
			solve(a == L, u, DBC)
			assign(u0, u)
			ar_max[i] = u0.vector().max()
			i=i+1
			file << (u, t)

		print ar_max.max()