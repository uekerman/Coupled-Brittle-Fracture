import numpy as np
import time
import precice


n = 20
dn = 1 / n

# generate mesh
y = np.linspace(0, 1, n + 1)

# preCICE setup
participant_name = "Dummy"
config_file_name = "precice-config.xml"
solver_process_index = 0
solver_process_size = 1
interface = precice.Interface(participant_name, config_file_name, solver_process_index, solver_process_size)

mesh_name = "Dummy-Mesh"
mesh_id = interface.get_mesh_id(mesh_name)

lmbda_id = interface.get_data_id("Lmbda", mesh_id)
gc_id = interface.get_data_id("Gc", mesh_id)

L = 1e-3
N = 100

axis = np.linspace(-L/2,L/2,N+1)

positions = np.transpose([np.tile(axis, len(axis)), np.repeat(axis, len(axis))])

vertex_ids = interface.set_mesh_vertices(mesh_id, positions)

precice_dt = interface.initialize()

step = 0

while interface.is_coupling_ongoing():
  
  print("Generating data")
  time.sleep(0.2)
  lmbda = np.ones((N+1)*(N+1)) * 121153.8e6 # First Lamé parameter
  gc = np.ones((N+1)*(N+1)) * 2.7e3 # Fracture toughness.
  for i in range(int(N*3/4), N+1):
    gc[(N+1)*int(N/2) + i] *= (1.0 - step / 100)
  
  interface.write_block_scalar_data(lmbda_id, vertex_ids, lmbda)
  interface.write_block_scalar_data(gc_id, vertex_ids, gc)
  
  precice_dt = interface.advance(precice_dt)
  
  step = step + 1
  
interface.finalize()

