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
        dh = 1.5 * mm
        dc = 0.35 * mm
        w0 = 55.49 * mm
        l0 = 42.99 * mm

        pcb_oversize = 2.0
        w1 = w0 * pcb_oversize
        l1 = l0 * pcb_oversize

        air_oversize = 2.0
        w2 = w1 * air_oversize
        l2 = l1 * air_oversize
        h2 = 50.0 * dh * air_oversize

        pml_oversize = 1.5
        w3 = w2 * pml_oversize
        l3 = l2 * pml_oversize
        h3 = h2 * pml_oversize

        w_feed = 2.0 * mm

        self.dims = {}
        self.dims['dh'] = dh
        self.dims['dc'] = dc
        self.dims['w0'] = w0
        self.dims['l0'] = l0
        self.dims['w1'] = w1
        self.dims['l1'] = l1
        self.dims['w_feed'] = w_feed
        self.dims['w2'] = w2
        self.dims['l2'] = l2
        self.dims['h2'] = h2
        self.dims['w3'] = w3
        self.dims['l3'] = l3
        self.dims['h3'] = h3

        self.tags = {}

        gmsh.initialize()
        gmsh.model.add(self.name)
        self._create_antenna()
        self._create_environment()
        gmsh.model.occ.synchronize()
        self._set_mesh_settings()
        self._create_groups()

        # gmsh.model.mesh.generate(3)
        # gmsh.write(self.name + '.msh')

    def _create_antenna(self):
        occ = gmsh.model.occ

        dh = self.dims['dh']
        dc = self.dims['dc']
        w0 = self.dims['w0']
        l0 = self.dims['l0']
        w1 = self.dims['w1']
        l1 = self.dims['l1']
        w_feed = self.dims['w_feed']

        # patch with feeder
        tag = occ.addRectangle(-w0 * 0.5, -l0 * 0.5, dh, w0, l0)
        patch2d = (2, tag)
        tag = occ.addRectangle(-w_feed * 0.5, l0 * 0.5,
                               dh, w_feed, (l1 - l0) * 0.5)
        feed2d = (2, tag)
        tags, _ = occ.fuse([patch2d], [feed2d], 0)
        patch2d = tags[0]
        tags = occ.extrude([patch2d], 0.0, 0.0, dc)
        patch3d = tags[1]

        # substrate rect
        tag = occ.addRectangle(-w1 * 0.5, -l1 * 0.5, 0.0, w1, l1)
        pcb2d = (2, tag)
        # port points on substrate rect
        tag = occ.addPoint(-w_feed * 0.5, l1 * 0.5, 0.0)
        port0d0 = (0, tag)
        tag = occ.addPoint(w_feed * 0.5, l1 * 0.5, 0.0)
        port0d1 = [0, tag]
        tags, _ = occ.fragment([pcb2d], [port0d0, port0d1])
        # extrude substrate and port
        pcb2d = tags[0]
        tags = occ.extrude([pcb2d], 0.0, 0.0, dh)
        pcb3d = tags[1]

        # extrude gnd foil
        tags = occ.extrude([pcb2d], 0.0, 0.0, -dc)
        gnd3d = tags[1]

        occ.synchronize()
        copper_skin = gmsh.model.getBoundary([gnd3d, patch3d])

        self.tags['pcb2d'] = pcb2d
        self.tags['pcb3d'] = pcb3d
        self.tags['patch2d'] = patch2d
        self.tags['patch3d'] = patch3d
        self.tags['gnd3d'] = gnd3d
        self.tags['copper_skin'] = copper_skin

    def _create_environment(self):
        occ = gmsh.model.occ

        dh = self.dims['dh']
        w2 = self.dims['w2']
        l2 = self.dims['l2']
        h2 = self.dims['h2']
        w3 = self.dims['w3']
        l3 = self.dims['l3']
        h3 = self.dims['h3']

        pcb2d = self.tags['pcb2d']
        pcb3d = self.tags['pcb3d']
        patch2d = self.tags['patch2d']
        patch3d = self.tags['patch3d']
        gnd3d = self.tags['gnd3d']

        # airbox
        tag = occ.addBox(-w2 * 0.5, -l2 * 0.5, (dh - h2) * 0.5, w2, l2, h2)
        air3d = (3, tag)
        tags, _ = occ.cut([air3d], [pcb3d, gnd3d, patch3d], 0, True, False)
        air3d = tags[0]
        # pml - perfect matched layer
        tag = occ.addBox(-w3 * 0.5, -l3 * 0.5, (dh - h3) * 0.5, w3, l3, h3)
        pml3d = (3, tag)
        occ.synchronize()
        pmlbox = gmsh.model.getBoundary([pml3d])

        tags, _ = occ.cut([pml3d], [air3d], 0, True, False)
        pml3d = tags[0]
        occ.remove([tags[1]], recursive=True)

        self.tags['air3d'] = air3d
        self.tags['pml3d'] = pml3d
        self.tags['pmlbox'] = pmlbox

    def _set_mesh_settings(self):
        # general settings
        gmsh.option.setNumber('Mesh.Algorithm3D', 4)  # (1=Delaunay, 4=Frontal)
        gmsh.option.setNumber('Mesh.Optimize', 1)
        gmsh.option.setNumber('Mesh.Smoothing', 5)

        # mesh sizes by elements
        mm = 1.0e-3
        mesh_size_condutor = 1.0 * mm
        mesh_size_substrate = 1.0 * mm
        mesh_size_environment = 10.0 * mm
        air3d = self.tags['air3d']
        gnd3d = self.tags['gnd3d']
        patch3d = self.tags['patch3d']
        pcb3d = self.tags['pcb3d']
        pml3d = self.tags['pml3d']

        tags = gmsh.model.getBoundary([air3d, pml3d], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_environment)

        tags = gmsh.model.getBoundary([pcb3d], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_substrate)

        tags = gmsh.model.getBoundary([patch3d, gnd3d], False, False, True)
        gmsh.model.mesh.setSize(tags, mesh_size_condutor)

    def _create_groups(self):
        dh = self.dims['dh']
        dc = self.dims['dc']
        l1 = self.dims['l1']
        w_feed = self.dims['w_feed']

        air3d = self.tags['air3d']
        gnd3d = self.tags['gnd3d']
        patch2d = self.tags['patch2d']
        patch3d = self.tags['patch3d']
        pcb2d = self.tags['pcb2d']
        pcb3d = self.tags['pcb3d']
        pml3d = self.tags['pml3d']
        pmlbox = self.tags['pmlbox']
        copper_skin = self.tags['copper_skin']

        eps = 1.0e-4
        tags = gmsh.model.occ.getEntitiesInBoundingBox(
            -w_feed * 0.5 - eps,
            l1 * 0.5 - eps,
            dh - eps,
            w_feed * 0.5 + eps,
            l1 * 0.5 + eps,
            dh + dc + eps,
            2
        )

        tag_feed = tags[0][1]
        tag = gmsh.model.addPhysicalGroup(2, [tag_feed])
        gmsh.model.setPhysicalName(2, tag, 'SkinFeed')

        ids = np.array(copper_skin)[:, 1]
        tag = gmsh.model.addPhysicalGroup(2, ids)
        gmsh.model.setPhysicalName(2, tag, 'SkinConductor')

        # tag = gmsh.model.getBoundary([pcb3d])
        # tag = np.array(tag)[:, 1]
        # tag = tag[tag != tag_feed]
        # tag = gmsh.model.addPhysicalGroup(2, tag)
        # gmsh.model.setPhysicalName(2, tag, 'Substrate')
        tag = gmsh.model.addPhysicalGroup(3, [pcb3d[1]])
        gmsh.model.setPhysicalName(3, tag, 'Substrate')

        tag = gmsh.model.addPhysicalGroup(3, [air3d[1]])
        gmsh.model.setPhysicalName(3, tag, 'Air')

        tag = gmsh.model.addPhysicalGroup(3, [pml3d[1]])
        gmsh.model.setPhysicalName(3, tag, 'Pml')

        ids = np.array(pmlbox)[:, 1]
        tag = gmsh.model.addPhysicalGroup(2, ids)
        gmsh.model.setPhysicalName(2, tag, 'SigmaInf')


# mspa = Mspa('mspa')
# if '-nopopup' not in sys.argv:
#     gmsh.fltk.run()

# gmsh.finalize()
