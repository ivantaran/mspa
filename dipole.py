
from gmsh import model
from gmsh import onelab
from gmsh import option
import gmsh
import numpy as np
import sys


occ = model.occ
field = model.mesh.field


class SimpleDipole(object):
    '''
    Simple dipole antenna
    '''

    def __init__(self, name='untitled'):
        super().__init__()
        self.name = name
        self.dims = {}
        self.tags = {}
        self._length = 0.488
        self._r = 0.01
        self._gap = 0.01
        self.refresh()

    @property
    def length(self):
        """Rod length"""
        return self._length

    @property
    def gap(self):
        """Gap"""
        return self._gap

    @property
    def r(self):
        """Rod radius"""
        return self._r

    @length.setter
    def length(self, value):
        self._length = value
        self.refresh()

    @gap.setter
    def gap(self, value):
        self._gap = value
        self.refresh()

    @r.setter
    def r(self, value):
        self._r = value
        self.refresh()

    def refresh(self):
        gmsh.clear()

        model.add(self.name)
        self._create_antenna()
        occ.synchronize()
        self._set_mesh_settings()
        self._create_groups()

    def _create_antenna(self):

        length = self.length
        gap = self.gap
        x = gap * 0.5
        r = self.r
        tag = occ.add_cylinder(-x, 0.0, 0.0, gap, 0.0, 0.0, r)
        vol_feed = (3, tag)
        tag = occ.add_cylinder(x, 0.0, 0.0, length, 0.0, 0.0, r)
        vol_cyl1 = (3, tag)
        tag = occ.add_cylinder(-x, 0.0, 0.0, -length, 0.0, 0.0, r)
        vol_cyl2 = (3, tag)

        air_r = 1.3
        pml_r = air_r + 0.2
        tag = occ.add_sphere(0.0, 0.0, 0.0, air_r)
        vol_air = (3, tag)
        tag = occ.add_sphere(0.0, 0.0, 0.0, pml_r)
        vol_pml = (3, tag)

        occ.synchronize()
        tags, _ = occ.cut([vol_pml], [vol_air],
                          tag=0, removeObject=True, removeTool=False)
        vol_pml = tags[0]

        occ.synchronize()
        tags, _ = occ.cut([vol_air], [vol_cyl1, vol_cyl2, vol_feed],
                          tag=0, removeObject=True, removeTool=False)
        vol_air = tags[0]

        occ.synchronize()
        occ.remove_all_duplicates()

        self.tags['vol_cond'] = [vol_cyl1, vol_cyl2]
        self.tags['vol_feed'] = vol_feed
        self.tags['vol_air'] = vol_air
        self.tags['vol_pml'] = vol_pml

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
        length = self.length + self.gap * 0.5
        r = self.r

        field.add("Distance", 1)
        field.set_numbers(1, "SurfacesList", [14, 15, 18])
        field.set_number(1, "NumPointsPerCurve", 100)

        field.add("Threshold", 2)
        field.set_number(2, "InField", 1)
        field.set_number(2, "SizeMin", r)
        field.set_number(2, "SizeMax", 0.2)
        field.set_number(2, "DistMin", 0.01)
        field.set_number(2, "DistMax", 0.01)

        field.set_as_background_mesh(2)

    def _create_groups(self):

        # sur_feed = self.tags['sur_feed']
        # sur_conductor = self.tags['sur_conductor']
        vol_feed = self.tags['vol_feed']
        vol_air = self.tags['vol_air']
        vol_pml = self.tags['vol_pml']
        vol_cond = self.tags['vol_cond']

        tags = model.get_boundary([vol_feed])
        sur_feed = tags[0]
        tag = model.add_physical_group(2, [sur_feed[1]])
        model.set_physical_name(2, tag, 'SkinFeed')
        model.set_color([vol_feed], 0, 0, 255, 128, True)

        tag = model.add_physical_group(3, [vol_air[1], vol_feed[1]])
        model.set_physical_name(3, tag, 'Air')
        model.set_color([vol_air], 255, 0, 0, 128, True)

        tag = model.add_physical_group(3, [vol_pml[1]])
        model.set_physical_name(3, tag, 'Pml')
        model.set_color([vol_pml], 255, 0, 0, 128, True)

        tags = model.get_boundary([vol_pml])
        sur_pml = tags[0]
        tag = model.add_physical_group(2, [sur_pml[1]])
        model.set_physical_name(2, tag, 'SigmaInf')

        sur_conductor = np.array([], dtype=int)
        for vol in vol_cond:
            tags = model.get_boundary([vol])
            tags = np.array(tags)
            sur_conductor = np.hstack((sur_conductor, tags[:, 1]))
        sur_conductor = np.unique(sur_conductor)
        tag = model.add_physical_group(2, sur_conductor)
        model.set_physical_name(2, tag, 'SkinConductor')
        model.set_color(vol_cond, 202, 204, 206, 128, True)


MODEL_NAME = 'dipole'
gmsh.initialize()
antenna = SimpleDipole(MODEL_NAME)
gmsh.open('fullwave.pro')
model.set_current(MODEL_NAME)
model.mesh.generate(3)
# gmsh.write(f'{MODEL_NAME}.msh')
gmsh.write('fullwave.msh')
onelab.run()

if "-nopopup" not in sys.argv:
    gmsh.fltk.initialize()
    gmsh.fltk.run()

gmsh.finalize()
