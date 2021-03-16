from fenics import *
from fenicsprecice import Adapter

mesh = Mesh('Gapmesh.xml')
x = mesh.coordinates()
x[:, :] *= 0.001 # mesh is scaled from m to mm to better fit with the existing fracture code

V = FunctionSpace(mesh, "Lagrange", 1)

def boundary(x):
    return x[0] < 0.001 + DOLFIN_EPS or x[0] > 0.0015 - DOLFIN_EPS
    
# define coupling domain as complete domain
def coupling(x):
    return True    

dim = 2
coupling_domain = AutoSubDomain(coupling)

# the FEniCS-preCICE adapter is configured via a json file. In there, we also need to specify some arbitrary (unused) read data.
# We are currently trying to get rid of this technical restriction.
precice = Adapter(adapter_config_filename="precice-adapter-config.json")

precice_dt = precice.initialize(coupling_domain, write_object=V) # here we define where the coupling should happen

    
# Define boundary condition on left and right boundary
u0 = Constant(0.0)
bc = DirichletBC(V, u0, boundary)

# Define variational problem, simply copied from FEniCS demo
u = TrialFunction(V)
v = TestFunction(V)
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)
a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

while precice.is_coupling_ongoing(): # potential time loop, currently only one timestep is done

  # Compute solution
  u = Function(V)
  solve(a == L, u, bc)
  
  # Save solution in VTK format
  file = File("poisson.pvd")
  file << u
  
  # write data to preCICE and advance coupling  
  precice.write_data(u)
  precice.advance(precice_dt)

precice.finalize()
