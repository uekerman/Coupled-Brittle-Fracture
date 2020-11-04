#!/usr/bin/env python3

from nutils import cli, types, mesh, function, solver, export, transform, topology
import numpy, numpy.random, typing, treelog, matplotlib.collections
import precice
from mpi4py import MPI


unit = types.unit(m=1, s=1, g=1e-3, N='kg*m/s2', Pa='N/m2')

def main(X:unit['m'], Y:unit['m'], l0:unit['m'], degree:int, du:unit['m']):

  '''
  Mechanical test case

  .. arguments::

     X [0.5mm]
       Domain size in x direction.

     Y [0.04mm]
       Domain size in y direction.

     l0 [0.015mm]
       Charateristic length scale.    

     degree [1]
       Polynomial degree of the approximation.

     du [0.01mm]
       Applied displacement.
  '''

  assert degree > 0

  # create the mesh
  topo, geom = mesh.rectilinear([numpy.linspace(0.001,0.001+X,31), numpy.linspace(0.001,0.001+Y,11)])
    
  # prepare the integration and post processing samples
  ipoints = topo.sample('gauss', 2*degree)
  bezier  = topo.sample('bezier', 2*degree)

  # initialize the namespace
  ns        = function.Namespace()
  ns.x      = geom
  ns.ubasis = topo.basis('th-spline', degree=degree).vector(topo.ndims)
  ns.dbasis = topo.basis('th-spline', degree=degree)
  ns.Hbasis = ipoints.basis()
  ns.u_i    = 'ubasis_ni ?solu_n'
  ns.d      = 'dbasis_n  ?sold_n'
  ns.H0     = 'Hbasis_n ?solH0_n'
  ns.l0     = l0
  ns.du     = du
  
  # volume coupling fields
  ns.Gc    = 'dbasis_n  ?gcdofs_n'
  ns.lmbda = 'dbasis_n  ?lmbdadofs_n'
  ns.mu    = 'dbasis_n  ?mudofs_n'

  # formulation
  ns.strain_ij = '( u_i,j + u_j,i ) / 2'
  ns.stress_ij = 'lmbda strain_kk δ_ij + 2 mu strain_ij'
  ns.psi       = 'stress_ij strain_ij / 2'
  ns.H         = function.max(ns.psi, ns.H0)
  ns.gamma     = '( d^2 + l0^2 d_,i d_,i ) / (2 l0)'

  # boundary condition for displacement field
  sqru  = topo.boundary['top'].integral('( u_i n_i - du )^2 d:x' @ ns, degree=degree*2)
  sqru += topo.boundary['bottom'].integral('( u_i n_i )^2 d:x' @ ns, degree=degree*2)
  sqru += topo.boundary['bottom'].boundary['left'].integral('u_i u_i d:x' @ ns, degree=degree*2)
  consu = solver.optimize('solu', sqru, droptol=1e-12)

  # initialize the solution vectors
  solu = numpy.zeros(ns.ubasis.shape[0])
  sold = numpy.zeros(ns.dbasis.shape[0])
  solH0 = ipoints.eval(0.)

  # preCICE setup
  configFileName = "precice-config.xml"
  participantName = "BrittleFracture"
  solverProcessIndex = 0
  solverProcessSize = 1
  interface = precice.Interface(participantName, configFileName, solverProcessIndex, solverProcessSize)

  # define coupling mesh
  meshName = "BrittleFracture-Mesh" 
  meshID = interface.get_mesh_id(meshName)
  couplingsample = topo.sample('gauss', degree=degree*2)
  vertices = couplingsample.eval(ns.x)
  dataIndices = interface.set_mesh_vertices(meshID, vertices)
  
  lmbda = 121153.8e6 # First Lamé parameter in Pa
  mu    = 80769.2e6 # Second Lamé parameter in Pa
  
  sqrl = couplingsample.integral((ns.lmbda - lmbda)**2)
  lmbdadofs = solver.optimize('lmbdadofs', sqrl, droptol=1e-12)
  
  sqrm = couplingsample.integral((ns.mu - mu)**2)
  mudofs = solver.optimize('mudofs', sqrm, droptol=1e-12)

  # coupling data
  gcID = interface.get_data_id("Gc", meshID)
  
  precice_dt = interface.initialize() # pseudo timestep size handled by preCICE

  nstep = 10000 # very high number of steps, end of simulation is steered by preCICE instead 
  
  # time loop
  with treelog.iter.fraction('step', range(nstep)) as counter:
    for istep in counter:
    
      if not interface.is_coupling_ongoing():
        break
      
      if interface.is_read_data_available():  
        gc = interface.read_block_scalar_data(gcID, dataIndices)
        gc_function = couplingsample.asfunction(gc)
        sqrg = couplingsample.integral((ns.Gc - gc_function)**2)
        gcdofs = solver.optimize('gcdofs', sqrg, droptol=1e-12)

      ############################
      # Phase field problem      #
      ############################

      resd  = ipoints.integral('( Gc / l0 ) ( d dbasis_n + l0^2 d_,i dbasis_n,i ) d:x' @ ns)
      resd += ipoints.integral('2 H ( d - 1 ) dbasis_n d:x' @ ns)

      sold = solver.solve_linear('sold', resd, arguments={'solu':solu, 'solH0':solH0, 'lmbdadofs':lmbdadofs, 'mudofs':mudofs, 'gcdofs':gcdofs})

      ############################
      # Elasticity problem       #
      ############################

      resu = topo.integral('( 1 - d )^2 ubasis_ni,j stress_ij d:x' @ ns, degree=2*degree)
      solu = solver.solve_linear('solu', resu, arguments={'sold':sold, 'lmbdadofs':lmbdadofs, 'mudofs':mudofs}, constrain=consu)

      # Update zero state and history field
      solH0 = ipoints.eval(ns.H, arguments={'solu':solu, 'solH0':solH0, 'lmbdadofs':lmbdadofs, 'mudofs':mudofs})
      
      # do the coupling
      precice_dt = interface.advance(precice_dt)

      ############################
      # Output                   #
      ############################
      
      # element-averaged history field
      transforms = ipoints.transforms[0]
      indicator  = function.kronecker(1., axis=0, length=len(transforms), pos=function.TransformsIndexWithTail(transforms, function.TRANS).index)
      areas, integrals = ipoints.integrate([indicator, indicator * ns.H], arguments={'solu':solu, 'solH0':solH0, 'lmbdadofs':lmbdadofs, 'mudofs':mudofs, 'gcdofs':gcdofs})
      H = indicator.dot(integrals/areas)

      # evaluate fields
      points, dvals, uvals, lvals, mvals, gcvals = bezier.eval(['x_i', 'd', 'u_i', 'lmbda', 'mu', 'Gc'] @ ns, arguments={'solu':solu, 'sold':sold, 'solH0':solH0, 'lmbdadofs':lmbdadofs, 'mudofs':mudofs, 'gcdofs':gcdofs})
      Hvals = bezier.eval(H, arguments={'solu':solu, 'solH0':solH0})
      
      with treelog.add(treelog.DataLog()):
        export.vtk('Solid_' + str(istep), bezier.tri, points, Gc=gcvals, D=dvals, U=uvals, H=Hvals)

  interface.finalize()

cli.run(main)
