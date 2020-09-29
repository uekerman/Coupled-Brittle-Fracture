import numpy as np
import time
import precice

# preCICE setup
participant_name = "Dummy"
config_file_name = "precice-config.xml"
solver_process_index = 0
solver_process_size = 1
interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

mesh_name = "Dummy-Mesh"
mesh_id   = interface.get_mesh_id(mesh_name)

lmbda_id = interface.get_data_id("Lmbda", mesh_id)
mu_id    = interface.get_data_id("Mu", mesh_id)
gc_id    = interface.get_data_id("Gc", mesh_id)


# define coupling mesh
L = 1e-3 # domain size in m
N = 100 # number of cells in each direction
axis = np.linspace(-L/2,L/2,N+1) # coordinates along one axis
positions = np.transpose([np.tile(axis, len(axis)), np.repeat(axis, len(axis))]) # mesh coordinates
vertex_ids = interface.set_mesh_vertices(mesh_id, positions)



precice_dt = interface.initialize() # pseudo timestep size handled by preCICE

step  = 0
nstep = 100

while interface.is_coupling_ongoing():
  
  print("Generating data")
  time.sleep(1.0) # for better readability of shell output
  lmbda = np.ones((N+1)*(N+1)) * 121153.8e6 # First Lamé parameter in Pa, constant everywhere
  mu    = np.ones((N+1)*(N+1)) * 80769.2e6 # Second Lamé parameter in Pa, constant everywhere
  gc    = np.ones((N+1)*(N+1)) * 2.7e3 # Fracture toughness in N/m, constant everywhere
  
  # modify fracture toughness at some arbitrary location (here vertical line, verically centered, right quarter of domain) 
  for i in range(int(N*3/4), N+1):
    gc[(N+1)*int(N/2) + i] *= (1.0 - step / nstep) # fracture toughness slowly decreases

  # write data to preCICE  
  interface.write_block_scalar_data(lmbda_id, vertex_ids, lmbda)
  interface.write_block_scalar_data(mu_id, vertex_ids, mu)
  interface.write_block_scalar_data(gc_id, vertex_ids, gc)
  
  # do the coupling
  precice_dt = interface.advance(precice_dt)
  
  step = step + 1
  
interface.finalize()

