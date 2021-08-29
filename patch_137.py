import gmsh
import math
import os
import sys
import numpy as np


class Mspa(object):
    '''
    Microstrip patch antenna
    '''

    def __init__(self, name='untitled'):
        super().__init__()
        self.name = name
        self.dims = {}
        self.tags = {}
        self.refresh()

    def refresh(self):
        gmsh.clear()
        mm = 1.0e-3

        r_cut = gmsh.onelab.get_number('Model/CutRadius')[0]
        size = gmsh.onelab.get_number('Model/PatchSize')[0]
        d_feed = gmsh.onelab.get_number('Model/FeedDistance')[0]
        # 200.0 * mm  # 240 * mm
        r_cut *= mm
        size *= mm
        d_feed *= mm
        d = 29.0 * mm
        w_path = size  # 640 791.0 * mm
        l_patch = size
        w_sub = 1170.0 * mm
        l_sub = 1170.0 * mm
        r_feed = 1.6 * mm
        r_shield = 4.0 * mm

        self.dims['d'] = d
        self.dims['gap'] = r_shield - r_feed
        self.dims['r_cut'] = r_cut
        self.dims['w_path'] = w_path
        self.dims['l_patch'] = l_patch
        self.dims['w_sub'] = w_sub
        self.dims['l_sub'] = l_sub
        self.dims['r_feed'] = r_feed
        self.dims['r_shield'] = r_shield
        self.dims['d_feed'] = d_feed

        gmsh.model.add(self.name)
        self._create_antenna()
        gmsh.model.occ.synchronize()
        self._set_mesh_settings()
        self._create_groups()

    def _create_antenna(self):
        occ = gmsh.model.occ

        d = self.dims['d']
        r_cut = self.dims['r_cut']
        w_path = self.dims['w_path']
        l_patch = self.dims['l_patch']
        w_sub = self.dims['w_sub']
        l_sub = self.dims['l_sub']
        r_feed = self.dims['r_feed']
        d_feed = self.dims['d_feed']
        r_shield = self.dims['r_shield']

        # substrate rect
        tag = occ.add_box(-0.5 * w_sub, -0.5 * l_sub, -
                          0.5 * d, w_sub, l_sub, d)
        vol_substrate = (3, tag)

        tag = occ.add_box(-0.5 * w_path, -0.5 * l_patch,
                          -0.5 * d, w_path, l_patch, d)
        vol_patch = (3, tag)

        tag = occ.add_cylinder(-0.5 * w_path, 0.5 *
                               l_patch, -0.5 * d, 0.0, 0.0, d, r_cut)
        vol_cyl1 = (3, tag)
        tag = occ.add_cylinder(0.5 * w_path, -0.5 *
                               l_patch, -0.5 * d, 0.0, 0.0, d, r_cut)
        vol_cyl2 = (3, tag)

        tag = occ.add_cylinder(0.0, 0.0, -0.5 * d, 0.0, 0.0, d, r_cut/2)
        vol_cyl3 = (3, tag)

        tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        vol_feed = (3, tag)
        tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        vol_feed = (3, tag)
        tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_shield)
        vol_shield = (3, tag)

        # w0 = 0.05
        # l0 = 0.04
        # l1 = 0.40
        # tag = occ.add_box(-0.5 * w_path, -0.5 * l0, -0.5 * d, w0, l0, d)
        # vol_1 = (3, tag)
        # tag = occ.add_box(-0.5 * w_path + w0, -0.5 * l1, -0.5 * d, l0, l1, d)
        # vol_2 = (3, tag)

        tags, _ = occ.cut([vol_patch], [vol_cyl1, vol_cyl2, vol_cyl3, vol_shield, vol_feed],
                          tag=0, removeObject=True, removeTool=True)
        vol_patch = tags[0]

        tag = occ.addSphere(0.0, 0.0, 0.0, l_sub)
        vol_air = (3, tag)
        tag = occ.addSphere(0.0, 0.0, 0.0, l_sub + 0.20)
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
        vol_substrate1 = tags[0]
        vol_substrate2 = tags[2]
        # vol_substrate2 = tags[1]

        occ.synchronize()
        occ.removeAllDuplicates()

        self.tags['sur_feed'] = (2, 44)
        self.tags['sur_conductor'] = [38, 28, 42, 24, 26, 43, 23, 22, 33]
        # self.tags['sur_feed'] = (2, 36)
        # self.tags['sur_conductor'] = [19, 20, 21, 25, 23, 34]

        self.tags['sur_pml'] = sur_pml
        self.tags['vol_air'] = vol_air
        self.tags['vol_pml'] = vol_pml
        self.tags['vol_substrate'] = [vol_patch[1],
                                      vol_substrate1[1], vol_substrate2[1], 14]

    def _set_mesh_settings(self):
        gmsh.option.setNumber('General.Antialiasing', 1)
        gmsh.option.setNumber('General.AlphaBlending', 1)
        gmsh.option.setNumber('View.FakeTransparency', 1)

        # 1: MeshAdapt, 2: Automatic, 3: Initial mesh only,
        # 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG,
        # 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
        gmsh.option.setNumber('Mesh.Algorithm', 2)
        # 1: Delaunay, 3: Initial mesh only,
        # 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
        gmsh.option.setNumber('Mesh.Algorithm3D', 4)
        gmsh.option.setNumber('Mesh.Optimize', 1)
        gmsh.option.setNumber('Mesh.Smoothing', 5)
        gmsh.option.setNumber('Mesh.SmoothNormals', 1)
        # 0: By element type
        # 1: By elementary entity
        # 2: By physical group
        # 3: By mesh partition
        gmsh.option.setNumber('Mesh.ColorCarousel', 2)
        gmsh.option.setNumber('Mesh.VolumeEdges', 0)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        # mesh sizes by elements
        d = self.dims['d']
        r_cut = self.dims['r_cut']
        w_path = self.dims['w_path']
        l_patch = self.dims['l_patch']
        w_sub = self.dims['w_sub']
        l_sub = self.dims['l_sub']
        r_feed = self.dims['r_feed']
        d_feed = self.dims['d_feed']
        r_shield = self.dims['r_shield']

        mm = 1.0e-3
        sur_feed = self.tags['sur_feed']
        sur_pml = self.tags['sur_pml']
        vol_air = self.tags['vol_air']
        vol_pml = self.tags['vol_pml']
        vol_substrate = self.tags['vol_substrate']

        tags = gmsh.model.getBoundary(
            [
                (3, vol_substrate[0]),
                (3, vol_substrate[1]),
                (3, vol_substrate[2])
            ], False, False, False)
        tags = gmsh.model.getBoundary(tags, False, False, False)
        a = np.array(tags)
        a = list(np.unique(a[:, 1]))
        # b = [43, 31, 33, 50, 32, 54]
        # c = [v for v in a if v not in b]

        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", a)
        gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 20)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", 0.01)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", 0.20)
        gmsh.model.mesh.field.setNumber(2, "DistMin", 0.00)
        gmsh.model.mesh.field.setNumber(2, "DistMax", 0.20)

        gmsh.model.mesh.field.add("Cylinder", 3)
        gmsh.model.mesh.field.setNumber(3, "Radius", 0.006)
        gmsh.model.mesh.field.setNumber(3, "VIn", 0.001)
        gmsh.model.mesh.field.setNumber(3, "VOut", 0.30)
        gmsh.model.mesh.field.setNumber(3, "XAxis", 0.00)
        gmsh.model.mesh.field.setNumber(3, "XCenter", 0.00)
        gmsh.model.mesh.field.setNumber(3, "YAxis", 0.00)
        gmsh.model.mesh.field.setNumber(3, "YCenter", -d_feed)
        gmsh.model.mesh.field.setNumber(3, "ZAxis", d)
        gmsh.model.mesh.field.setNumber(3, "ZCenter", 0.00)

        gmsh.model.mesh.field.add("Min", 4)
        gmsh.model.mesh.field.setNumbers(4, "FieldsList", [2, 3])

        gmsh.model.mesh.field.setAsBackgroundMesh(4)

    def _create_groups(self):

        sur_feed = self.tags['sur_feed']
        sur_conductor = self.tags['sur_conductor']
        sur_pml = self.tags['sur_pml']
        vol_air = self.tags['vol_air']
        vol_pml = self.tags['vol_pml']
        vol_substrate = self.tags['vol_substrate']

        tag = gmsh.model.addPhysicalGroup(2, [sur_feed[1]])
        gmsh.model.setPhysicalName(2, tag, 'SkinFeed')

        tag = gmsh.model.addPhysicalGroup(2, sur_conductor)
        gmsh.model.setPhysicalName(2, tag, 'SkinConductor')

        tag = gmsh.model.addPhysicalGroup(3, vol_substrate)
        gmsh.model.setPhysicalName(3, tag, 'Substrate')

        tag = gmsh.model.addPhysicalGroup(3, [vol_air[1]])
        gmsh.model.setPhysicalName(3, tag, 'Air')

        tag = gmsh.model.addPhysicalGroup(3, [vol_pml[1]])
        gmsh.model.setPhysicalName(3, tag, 'Pml')
        gmsh.model.setColor([vol_pml], 255, 0, 0, 16, True)

        tag = gmsh.model.addPhysicalGroup(2, [sur_pml[1]])
        gmsh.model.setPhysicalName(2, tag, 'SigmaInf')
