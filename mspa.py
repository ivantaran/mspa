import gmsh
import math
import os
import sys
import numpy as np


MODEL_NAME = 'mspa'

mm = 1.0e-3
dh = 1.5 * mm
dc = 0.35 * mm
w0 = 55.49 * mm
l0 = 42.99 * mm

pcb_oversize = 1.5
w1 = w0 * pcb_oversize
l1 = l0 * pcb_oversize

air_oversize = 1.5
w2 = w1 * air_oversize
l2 = l1 * air_oversize
h2 = 10.0 * dh * air_oversize

pml_oversize = 1.5
w3 = w2 * pml_oversize
l3 = l2 * pml_oversize
h3 = h2 * pml_oversize

w_feed = 2.0 * mm

mesh_size_condutor = 1.5 * mm
mesh_size_substrate = 1.5 * mm
mesh_size_environment = 10 * mm

gmsh.initialize()

# gmsh.open('mspa.pro')
# gmsh.onelab.run()

# gmsh.model.setCurrent('mspa.geo')
# path = os.path.dirname(os.path.abspath(__file__))
# gmsh.merge(os.path.join(path, 'mspa.pro'))

gmsh.model.add(MODEL_NAME + '.geo')
occ = gmsh.model.occ

patch2d = occ.addRectangle(-w0 * 0.5, -l0 * 0.5, dh, w0, l0)
feed2d = occ.addRectangle(-w_feed * 0.5, l0 * 0.5, dh, w_feed, (l1 - l0) * 0.5)
tags, _ = occ.fuse([(2, patch2d)], [(2, feed2d)], 0)
patch2d = tags[0][1]

tags = occ.extrude([(2, patch2d)], 0.0, 0.0, dc)
patch3d = tags[1][1]

pcb2d = occ.addRectangle(-w1 * 0.5, -l1 * 0.5, 0.0, w1, l1)

port_point0 = occ.addPoint(-w_feed * 0.5, l1 * 0.5, 0.0)
port_point1 = occ.addPoint(w_feed * 0.5, l1 * 0.5, 0.0)
tags, _ = occ.fragment([(2, pcb2d)], [(0, port_point0),
                                      (0, port_point1)])
pcb2d = tags[0][1]
tags = occ.extrude([(2, pcb2d)], 0.0, 0.0, dh)
pcb3d = tags[1][1]

tags = occ.extrude([(2, pcb2d)], 0.0, 0.0, -dc)
gnd3d = tags[1][1]

air3d = occ.addBox(-w2 * 0.5, -l2 * 0.5, (dh - h2) * 0.5, w2, l2, h2)
tags, _ = occ.cut([(3, air3d)], [(3, pcb3d), (3, gnd3d),
                                 (3, patch3d)], 0, True, False)
air3d = tags[0][1]

pml3d = occ.addBox(-w3 * 0.5, -l3 * 0.5, (dh - h3) * 0.5, w3, l3, h3)
tags, _ = occ.cut([(3, pml3d)], [(3, air3d)], 0, True, False)
pml3d = tags[0][1]
occ.remove([tags[1]], recursive=True)

occ.synchronize()

gmsh.option.setNumber('Mesh.Algorithm3D', 4)  # (1=Delaunay, 4=Frontal)
gmsh.option.setNumber('Mesh.Optimize', 1)
gmsh.option.setNumber('Mesh.Smoothing', 5)

tag = gmsh.model.addPhysicalGroup(2, [patch2d, pcb2d])
gmsh.model.setPhysicalName(2, tag, 'Conductor')

eps = 1.0e-4
tags = gmsh.model.occ.getEntitiesInBoundingBox(
    -w_feed * 0.5 - eps, l1 * 0.5 - eps, -eps, w_feed * 0.5 + eps, l1 * 0.5 + eps, dh + eps, 2)
tag = gmsh.model.addPhysicalGroup(2, [tags[0][1]])
gmsh.model.setPhysicalName(2, tag, 'SkinFeed')

tag = gmsh.model.addPhysicalGroup(3, [patch3d, gnd3d])
gmsh.model.setPhysicalName(3, tag, 'ConductorVolume')

tag = gmsh.model.addPhysicalGroup(3, [pcb3d])
gmsh.model.setPhysicalName(3, tag, 'Substrate')

tag = gmsh.model.addPhysicalGroup(3, [air3d])
gmsh.model.setPhysicalName(3, tag, 'Air')

tag = gmsh.model.addPhysicalGroup(3, [pml3d])
gmsh.model.setPhysicalName(3, tag, 'Pml')


tags = gmsh.model.getBoundary([(3, pcb3d)], False, False, True)
gmsh.model.mesh.setSize(tags, mesh_size_substrate)

tags = gmsh.model.getBoundary([(3, air3d), (3, pml3d)], False, False, True)
gmsh.model.mesh.setSize(tags, mesh_size_environment)

tags = gmsh.model.getBoundary(
    [(3, patch3d), (3, gnd3d)], False, False, True)
gmsh.model.mesh.setSize(tags, mesh_size_condutor)

gmsh.model.mesh.generate(3)
gmsh.write(MODEL_NAME + '.msh')

# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()

# gmsh.finalize()
