
from gmsh import model
from gmsh import onelab
from gmsh import option
import gmsh
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
        self._d_feed = 0.12
        self._r_cut = 0.001
        self._patch_size = 0.850
        self.refresh()

    @property
    def patch_size(self):
        """Patch size"""
        return self._patch_size

    @property
    def d_feed(self):
        """Feed distance"""
        return self._d_feed

    @property
    def r_cut(self):
        """Cut radius"""
        return self._r_cut

    @patch_size.setter
    def patch_size(self, value):
        self._patch_size = value
        self.refresh()

    @d_feed.setter
    def d_feed(self, value):
        self._d_feed = value
        self.refresh()

    @r_cut.setter
    def r_cut(self, value):
        self._r_cut = value
        self.refresh()

    def refresh(self):
        gmsh.clear()

        mm = 1.0e-3
        d = 29.0 * mm
        w_path = self.patch_size
        l_patch = self.patch_size
        w_sub = 1170.0 * mm
        l_sub = 1170.0 * mm
        r_feed = 0.40 * mm
        r_shield = 1.1 * mm

        self.dims['d'] = d
        self.dims['gap'] = r_shield - r_feed
        self.dims['r_cut'] = self.r_cut
        self.dims['w_path'] = w_path
        self.dims['l_patch'] = l_patch
        self.dims['w_sub'] = w_sub
        self.dims['l_sub'] = l_sub
        self.dims['r_feed'] = r_feed
        self.dims['r_shield'] = r_shield
        self.dims['d_feed'] = self.d_feed

        model.add(self.name)
        self._create_antenna()
        model.occ.synchronize()
        self._set_mesh_settings()
        self._create_groups()

    def _create_antenna(self):
        occ = model.occ

        d = self.dims['d']
        r_cut = self.dims['r_cut']
        w_path = self.dims['w_path']
        l_patch = self.dims['l_patch']
        w_sub = self.dims['w_sub']
        l_sub = self.dims['l_sub']
        r_feed = self.dims['r_feed']
        d_feed = self.d_feed
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

        # tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        # vol_feed = (3, tag)
        # tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        # vol_feed = (3, tag)
        # tag = occ.add_cylinder(0.0, -d_feed, -0.5 * d, 0.0, 0.0, d, r_shield)
        # vol_shield = (3, tag)

        tag = occ.add_cylinder(-d_feed, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        vol_feed = (3, tag)
        tag = occ.add_cylinder(-d_feed, -d_feed, -0.5 * d, 0.0, 0.0, d, r_feed)
        vol_feed = (3, tag)
        tag = occ.add_cylinder(-d_feed, -d_feed, -0.5 *
                               d, 0.0, 0.0, d, r_shield)
        vol_shield = (3, tag)

        w0 = 0.05  # 0.05
        l0 = 0.03  # 0.03
        l1 = 0.45  # 0.45
        tag = occ.add_box(-0.5 * w_path, -0.5 * l0, -0.5 * d, w0, l0, d)
        vol_1 = (3, tag)
        tag = occ.add_box(
            -0.5 * w_path + w0,
            -0.5 * l1,
            -0.5 * d, l0,
            l1, d)
        vol_2 = (3, tag)
        tag = occ.add_box(0.5 * w_path, -0.5 * l0, -0.5 * d, -w0, l0, d)
        vol_3 = (3, tag)
        tag = occ.add_box(
            0.5 * w_path - w0,
            -0.5 * l1,
            -0.5 * d, -l0,
            l1, d)
        vol_4 = (3, tag)

        tags, _ = occ.cut([vol_patch], [vol_cyl1, vol_cyl2, vol_shield, vol_feed, vol_1, vol_2, vol_3, vol_4],
                          tag=0, removeObject=True, removeTool=True)
        vol_patch = tags[0]

        tag = occ.add_sphere(0.0, 0.0, 0.0, l_sub)
        vol_air = (3, tag)
        tag = occ.add_sphere(0.0, 0.0, 0.0, l_sub + 0.20)
        vol_pml = (3, tag)

        occ.synchronize()
        tags = model.get_boundary([vol_pml])
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
        vol_substrate2 = tags[1]

        occ.synchronize()
        occ.remove_all_duplicates()

        self.tags['sur_feed'] = (2, 52)
        self.tags['sur_conductor'] = [19, 20, 21, 23, 25, 50, 51]

        self.tags['sur_pml'] = sur_pml
        self.tags['vol_air'] = vol_air
        self.tags['vol_pml'] = vol_pml
        self.tags['vol_substrate'] = [vol_patch[1],
                                      vol_substrate1[1], vol_substrate2[1]]

    def _set_mesh_settings(self):
        option.set_number('General.Antialiasing', 1)
        option.set_number('General.AlphaBlending', 1)
        option.set_number('View.FakeTransparency', 1)

        # 1: MeshAdapt, 2: Automatic, 3: Initial mesh only,
        # 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG,
        # 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms
        option.set_number('Mesh.Algorithm', 2)
        # 1: Delaunay, 3: Initial mesh only,
        # 4: Frontal, 7: MMG3D, 9: R-tree, 10: HXT
        option.set_number('Mesh.Algorithm3D', 4)
        option.set_number('Mesh.Optimize', 1)
        option.set_number('Mesh.Smoothing', 5)
        option.set_number('Mesh.SmoothNormals', 1)
        # 0: By element type
        # 1: By elementary entity
        # 2: By physical group
        # 3: By mesh partition
        option.set_number('Mesh.ColorCarousel', 2)
        option.set_number('Mesh.VolumeEdges', 0)

        option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)
        option.set_number("Mesh.MeshSizeFromPoints", 0)
        option.set_number("Mesh.MeshSizeFromCurvature", 0)

        # mesh sizes by elements
        d = self.dims['d']
        d_feed = self.dims['d_feed']

        vol_substrate = self.tags['vol_substrate']

        tags = model.get_boundary(
            [
                (3, vol_substrate[0]),
                (3, vol_substrate[1]),
                (3, vol_substrate[2])
            ], False, False, False)

        tags = model.get_boundary(tags, False, False, False)
        a = np.array(tags)
        a = list(np.unique(a[:, 1]))

        model.mesh.field.add("Distance", 1)
        model.mesh.field.set_numbers(1, "CurvesList", a)
        model.mesh.field.set_number(1, "NumPointsPerCurve", 20)

        model.mesh.field.add("Threshold", 2)
        model.mesh.field.set_number(2, "InField", 1)
        model.mesh.field.set_number(2, "SizeMin", 0.01)
        model.mesh.field.set_number(2, "SizeMax", 0.20)
        model.mesh.field.set_number(2, "DistMin", 0.00)
        model.mesh.field.set_number(2, "DistMax", 0.20)

        model.mesh.field.add("Cylinder", 3)
        model.mesh.field.set_number(3, "Radius", 0.0011)
        model.mesh.field.set_number(3, "VIn", 0.0005)
        model.mesh.field.set_number(3, "VOut", 0.30)
        model.mesh.field.set_number(3, "XAxis", 0.00)
        model.mesh.field.set_number(3, "XCenter", -d_feed)
        # model.mesh.field.set_number(3, "XCenter", 0.00)
        model.mesh.field.set_number(3, "YAxis", 0.00)
        model.mesh.field.set_number(3, "YCenter", -d_feed)
        model.mesh.field.set_number(3, "ZAxis", d)
        model.mesh.field.set_number(3, "ZCenter", 0.00)

        model.mesh.field.add("Min", 4)
        model.mesh.field.set_numbers(4, "FieldsList", [2, 3])

        model.mesh.field.set_as_background_mesh(4)

    def _create_groups(self):

        sur_feed = self.tags['sur_feed']
        sur_conductor = self.tags['sur_conductor']
        sur_pml = self.tags['sur_pml']
        vol_air = self.tags['vol_air']
        vol_pml = self.tags['vol_pml']
        vol_substrate = self.tags['vol_substrate']

        tag = model.add_physical_group(2, [sur_feed[1]])
        model.set_physical_name(2, tag, 'SkinFeed')

        tag = model.add_physical_group(2, sur_conductor)
        model.set_physical_name(2, tag, 'SkinConductor')

        tag = model.add_physical_group(3, vol_substrate)
        model.set_physical_name(3, tag, 'Substrate')

        tag = model.add_physical_group(3, [vol_air[1]])
        model.set_physical_name(3, tag, 'Air')

        tag = model.add_physical_group(3, [vol_pml[1]])
        model.set_physical_name(3, tag, 'Pml')
        model.set_color([vol_pml], 255, 0, 0, 16, True)

        tag = model.add_physical_group(2, [sur_pml[1]])
        model.set_physical_name(2, tag, 'SigmaInf')
