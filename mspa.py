import gmsh
import math
import os
import sys
import numpy as np

# gmsh.open('mspa.pro')
# gmsh.onelab.run()

# gmsh.model.setCurrent('mspa.geo')
# path = os.path.dirname(os.path.abspath(__file__))
# gmsh.merge(os.path.join(path, 'mspa.pro'))

#


class Mspa(object):
    '''
    Microstrip patch antenna
    '''

    def __init__(self, name='untitled'):
        super().__init__()

        self.name = name

        mm = 1.0e-3
        mil = 0.0254 * mm

        d = 60.0 * mil
        w_line = 3.2 * mm
        w_path = 53.0 * mm
        l_patch = 52.0 * mm
        w_stub = 7.0 * mm
        l_stub = 15.5 * mm
        w_sub = 100.0 * mm
        l_sub = 100.0 * mm

        self.dims = {}
        self.dims['d'] = d
        self.dims['w_line'] = w_line
        self.dims['w_path'] = w_path
        self.dims['l_patch'] = l_patch
        self.dims['w_stub'] = w_stub
        self.dims['l_stub'] = l_stub
        self.dims['w_sub'] = w_sub
        self.dims['l_sub'] = l_sub

        self.tags = {}

        gmsh.initialize()
        gmsh.model.add(self.name)
        self._create_antenna()
        gmsh.model.occ.synchronize()
        self._set_mesh_settings()
        self._create_groups()

    def _create_antenna(self):
        occ = gmsh.model.occ

        d = self.dims['d']
        w_line = self.dims['w_line']
        w_path = self.dims['w_path']
        l_patch = self.dims['l_patch']
        w_stub = self.dims['w_stub']
        l_stub = self.dims['l_stub']
        w_sub = self.dims['w_sub']
        l_sub = self.dims['l_sub']

        # substrate rect
        tag = occ.addBox(-0.5 * w_sub, -0.5 * l_sub, -0.5 * d, w_sub, l_sub, d)
        vol_substrate = (3, tag)

        tag = occ.addRectangle(-0.5 * w_path, -0.5 *
                               l_patch, 0.5 * d, w_path, l_patch)
        sur_patch = (2, tag)

        tag = occ.addRectangle(-0.5 * w_stub, -0.5 *
                               l_stub, 0.5 * d, w_stub, l_stub)
        sur_stub1 = (2, tag)
        occ.translate([sur_stub1], 0.5 * (w_stub + w_line),
                      0.5 * (l_stub - l_patch), 0.0)

        tags = occ.copy([sur_stub1])
        sur_stub2 = tags[0]
        occ.translate([sur_stub2], -(w_stub + w_line), 0.0, 0.0)

        tags, _ = occ.cut([sur_patch], [sur_stub1, sur_stub2],
                          tag=0, removeObject=True, removeTool=True)
        sur_patch = tags[0]

        tags = occ.extrude([(1, 20)], 0.0, 0.0, -d)
        sur_feed = tags[1]

        tag = occ.addSphere(0.0, 0.0, 0.0, l_sub)
        vol_air = (3, tag)
        tag = occ.addSphere(0.0, 0.0, 0.0, l_sub * 1.2)
        vol_pml = (3, tag)

        occ.synchronize()
        tags = gmsh.model.getBoundary([vol_pml])
        sur_pml = tags[0]

        tags, _ = occ.cut([vol_pml], [vol_air],
                          tag=0, removeObject=True, removeTool=False)
        vol_pml = tags[0]

        tags, _ = occ.cut([vol_air], [vol_substrate],
                          tag=0, removeObject=True, removeTool=False)
        vol_air = tags[0]

        tags, _ = occ.cut([vol_substrate], [vol_patch],
                          tag=0, removeObject=True, removeTool=False)
        vol_substrate = tags[0]

        occ.synchronize()

        self.tags['sur_feed'] = (2, 18)
        self.tags['sur_gnd'] = (2, 24)
        self.tags['sur_gnd1'] = (2, 11)
        self.tags['sur_patch'] = (2, 9)
        self.tags['sur_pml'] = sur_pml
        self.tags['vol_air'] = vol_air
        self.tags['vol_patch'] = vol_patch
        self.tags['vol_pml'] = vol_pml
        self.tags['vol_substrate'] = vol_substrate

    def _set_mesh_settings(self):
        # general settings
        gmsh.option.setNumber('Mesh.Algorithm3D', 4)  # (1=Delaunay, 4=Frontal)
        gmsh.option.setNumber('Mesh.Optimize', 1)
        gmsh.option.setNumber('Mesh.Smoothing', 5)

        # mesh sizes by elements
        mm = 1.0e-3
        mesh_size_condutor = 1.5 * mm
        mesh_size_substrate = 1.5 * mm
        mesh_size_environment = 10.0 * mm
        sur_feed = self.tags['sur_feed']
        sur_gnd = self.tags['sur_gnd']
        sur_patch = self.tags['sur_patch']
        sur_pml = self.tags['sur_pml']
        vol_air = self.tags['vol_air']
        vol_patch = self.tags['vol_patch']
        vol_pml = self.tags['vol_pml']
        vol_substrate = self.tags['vol_substrate']

        tags = gmsh.model.getBoundary([vol_air, vol_pml], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_environment)

        tags = gmsh.model.getBoundary(
            [vol_substrate, vol_patch], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_substrate)

        tags = gmsh.model.getBoundary(
            [sur_feed, sur_gnd, sur_patch], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_condutor)

    def _create_groups(self):

        sur_feed = self.tags['sur_feed']
        sur_gnd = self.tags['sur_gnd']
        sur_gnd1 = self.tags['sur_gnd1']
        sur_patch = self.tags['sur_patch']
        sur_pml = self.tags['sur_pml']
        vol_air = self.tags['vol_air']
        vol_patch = self.tags['vol_patch']
        vol_pml = self.tags['vol_pml']
        vol_substrate = self.tags['vol_substrate']

        tag = gmsh.model.addPhysicalGroup(2, [sur_feed[1]])
        gmsh.model.setPhysicalName(2, tag, 'SkinFeed')

        tag = gmsh.model.addPhysicalGroup(
            2, [sur_gnd[1], sur_gnd1[1], sur_patch[1]])
        gmsh.model.setPhysicalName(2, tag, 'SkinConductor')

        tag = gmsh.model.addPhysicalGroup(3, [vol_substrate[1], vol_patch[1]])
        gmsh.model.setPhysicalName(3, tag, 'Substrate')

        tag = gmsh.model.addPhysicalGroup(3, [vol_air[1]])
        gmsh.model.setPhysicalName(3, tag, 'Air')

        tag = gmsh.model.addPhysicalGroup(3, [vol_pml[1]])
        gmsh.model.setPhysicalName(3, tag, 'Pml')

        tag = gmsh.model.addPhysicalGroup(2, [sur_pml[1]])
        gmsh.model.setPhysicalName(2, tag, 'SigmaInf')
