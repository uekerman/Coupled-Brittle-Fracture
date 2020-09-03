#!/usr/bin/env python3

from nutils import cli, types, mesh, function, solver, export, transform, topology
import numpy, numpy.random, typing, treelog, matplotlib.collections

unit = types.unit(m=1, s=1, g=1e-3, N='kg*m/s2', Pa='N/m2')

def main(L:unit['m'], c:unit['m'], l0:unit['m'], h0:unit['m'], nr:int, degree:int, ry0:unit['m'], lmbda:unit['Pa'], mu:unit['Pa'], Gc:unit['N/m'], nstep:int, du:unit['m']):

  '''
  Mechanical test case

  .. arguments::

     L [1.0mm]
       Domain size.

     c [0.5mm]
       Initial crack length.

     l0 [0.015mm]
       Charateristic length scale.

     h0 [0.25mm]
       Coarse-scale mesh size.     

     nr [4]
       Number of local refinements.

     degree [1]
       Polynomial degree of the approximation.

     ry0 [0.09mm]
       Half-size of the y-refinement zone.

     lmbda [121153.8MPa]
       First Lame parameter.

     mu [80769.2MPa]
       Second Lame parameter.

     Gc [2.7N/mm]
       Fracture toughness.

     nstep [150]
       Number of steps.

     du [0.0001mm]
       Increment in displacement.
  '''

  assert degree > 0

  m0 = int(L/h0) if int(L/h0)%2==0 else int(L/h0)+1
  m1 = (2**nr)*m0
  h1 = L/m1

  treelog.info('nelems (coarse)= {}'.format(m0))
  treelog.info('nelems (fine)  = {}'.format(m1))

  # create the mesh
  topo, geom = mesh.rectilinear([numpy.linspace(-L/2,L/2,m0+1)]*2)
  for ir in range(nr):
    bezier = topo.sample('vertex',0)
    vals = bezier.eval(geom,  separate=True)
    mask = [(abs(vals[i,1])<ry0).any() for i in bezier.index]
    topo = topo.refined_by(numpy.arange(len(topo),dtype=int)[mask])

  # create the initial fracture topology
  interfaces = topo.interfaces
  sbezier = interfaces.sample('uniform',1)
  vals = sbezier.eval(geom, separate=True)
  mask = numpy.concatenate([(abs(vals[i,1])<0.5*h1 and vals[i,0]<0) for i in sbezier.index])

  ctopo = topology.Topology(interfaces.references[mask], transforms=interfaces.transforms[mask], opposites=interfaces.opposites[mask])

  treelog.info('initial fracture length check: {:3.2f})'.format(ctopo.integrate(function.J(geom),ischeme='gauss1')))

  # prepare the integration and post processing samples
  ipoints = topo.sample('gauss', 2*degree)
  bezier  = topo.sample('bezier', 4)

  # initialize the namespace
  ns        = function.Namespace()
  ns.x      = geom
  ns.ubasis = topo.basis('th-spline', degree=degree).vector(topo.ndims)
  ns.dbasis = topo.basis('th-spline', degree=degree)
  ns.Hbasis = ipoints.basis()
  ns.u_i    = 'ubasis_ni ?solu_n'
  ns.d      = 'dbasis_n  ?sold_n'
  ns.H0     = 'Hbasis_n ?solH0_n'
  ns.lmbda  = lmbda
  ns.mu     = mu
  ns.l0     = l0
  ns.Gc     = Gc
  ns.du     = du

  # formulation
  ns.strain_ij = '( u_i,j + u_j,i ) / 2'
  ns.stress_ij = 'lmbda strain_kk δ_ij + 2 mu strain_ij'
  ns.psi       = 'stress_ij strain_ij / 2'
  ns.H         = function.max(ns.psi, ns.H0)
  ns.gamma     = '( d^2 + l0^2 d_,i d_,i ) / (2 l0)'

  # incremental solution constraints
  sqru  = topo.boundary['top'].integral('( u_i n_i - du )^2 d:x' @ ns, degree=degree*2)
  sqru += topo.boundary['bottom'].integral('( u_i n_i )^2 d:x' @ ns, degree=degree*2)
  sqru += topo.boundary['bottom'].boundary['left'].integral('u_i u_i d:x' @ ns, degree=degree*2)
  consu = solver.optimize('solu', sqru, droptol=1e-12)

  sqrd  = ctopo.integral('(d - 1)^2 d:x' @ ns, degree=degree*2)
  consd = solver.optimize('sold', sqrd, droptol=1e-12)

  # initialize the solution vectors
  solu = numpy.zeros(ns.ubasis.shape[0])
  sold = consd.copy()
  sold[numpy.isnan(consd)] = 0.
  solH0 = ipoints.eval(0.)

  # Initialize data containers
  data_clen  = []
  data_disp  = []
  data_force = []

  # Make loading steps
  with treelog.iter.fraction('step', range(nstep)) as counter:
    for istep in counter:

      ############################
      # Phase field problem      #
      ############################

      resd  = ipoints.integral('( Gc / l0 ) ( d dbasis_n + l0^2 d_,i dbasis_n,i ) d:x' @ ns)
      resd += ipoints.integral('2 H ( d - 1 ) dbasis_n d:x' @ ns)

      sold = solver.solve_linear('sold', resd, arguments={'solu':solu, 'solH0':solH0}, constrain=consd)

      ############################
      # Elasticity problem       #
      ############################

      resu = topo.integral('( 1 - d )^2 ubasis_ni,j stress_ij d:x' @ ns, degree=2*degree)
      solu = solver.solve_linear('solu', resu, arguments={'sold':sold}, constrain=istep*consu)

      # Update zero state and history field
      solH0 = ipoints.eval(ns.H, arguments={'solu':solu, 'solH0':solH0})

      ############################
      # Output                   #
      ############################

      clen = ipoints.integrate('gamma d:x' @ ns, arguments={'sold':sold})

      # append data containers
      blen, bdisp, bforc = topo.boundary['top'].integrate(['d:x', 'u_i n_i d:x', '( 1 - d )^2 stress_ij n_i n_j d:x'] @ ns(solu=solu, sold=sold), degree=2*degree)
      bdisp /= blen
      data_disp.append(bdisp/unit('1.0mm'))
      data_force.append(bforc*unit('1.0mm'))
      data_clen.append(clen/unit('1.0mm'))

      # element-averaged history field
      transforms = ipoints.transforms[0]
      indicator  = function.kronecker(1., axis=0, length=len(transforms), pos=function.TransformsIndexWithTail(transforms, function.TRANS).index)
      areas, integrals = ipoints.integrate([indicator, indicator * ns.H], arguments={'solu':solu, 'solH0':solH0})
      H = indicator.dot(integrals/areas)

      # evaluate fields
      points, dvals, uvals, psivals = bezier.eval(['x_i', 'd', 'u_i', 'psi'] @ ns, arguments={'solu':solu, 'sold':sold, 'solH0':solH0})
      Hvals = bezier.eval(H, arguments={'solu':solu, 'solH0':solH0})

      with export.mplfigure('solution.png') as fig:
        ax = fig.add_subplot(111)
        ax.plot(data_disp, data_force, '.-')
        ax.set_title('c = {:.2E}'.format(data_clen[-1]))
        ax.set_xlabel('Displacement')
        ax.set_ylabel('Force')
        ax.grid()

      with export.mplfigure('displacement.png') as fig:
        ax = fig.add_subplot(111)
        im = ax.tripcolor(points[:,0], points[:,1], bezier.tri, numpy.linalg.norm(uvals,axis=1), shading='gouraud', cmap='jet', rasterized=True)
        ax.add_collection(matplotlib.collections.LineCollection(points[bezier.hull], colors='k', linewidth=0.5, alpha=0.2))
        ax.autoscale(enable=True, axis='both', tight=False)
        fig.colorbar(im)

      with export.mplfigure('damage.png') as fig:
        ax = fig.add_subplot(111)
        im = ax.tripcolor(points[:,0], points[:,1], bezier.tri, dvals, shading='gouraud', cmap='jet', rasterized=True)
        ax.add_collection(matplotlib.collections.LineCollection(points[bezier.hull], colors='k', linewidth=0.5, alpha=0.2))
        ax.autoscale(enable=True, axis='both', tight=False)
        fig.colorbar(im)
        im.set_clim(0,1)

      with export.mplfigure('history.png') as fig:
        ax = fig.add_subplot(111)
        im = ax.tripcolor(points[:,0], points[:,1], bezier.tri, Hvals, shading='gouraud', cmap='jet', rasterized=True)
        ax.add_collection(matplotlib.collections.LineCollection(points[bezier.hull], colors='k', linewidth=0.5, alpha=0.2))
        ax.autoscale(enable=True, axis='both', tight=False)
        fig.colorbar(im)

cli.run(main)
