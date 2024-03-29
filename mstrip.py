
# from comsol_patch_1575 import Mspa
from gmsh import model
from gmsh import onelab
from gmsh import option
from gmsh import plugin
from matplotlib import pyplot
from numpy.lib.type_check import imag
from patch_137 import Mspa
from pprint import pprint
from pygetdp import Group, Function, Problem
from pygetdp.helpers import build_example_png, print_html
from scipy.constants import mu_0, epsilon_0, pi, speed_of_light
import gmsh
import numpy as np
import os
import sys


# GDICT1 = {
#     'Point': 1,
#     'Line': 2,
#     'Triangle': 3,
#     'Quadrangle': 4,
#     'Tetrahedron': 4,
#     'Hexahedron': 8,
#     'Prism': 6,
# }

# GDICT2 = GDICT1

GDICT1 = {
    'Point': 1,
    'Line': 3,
    'Triangle': 4,
    'Quadrangle': 4,
    'Tetrahedron': 4,
    'Hexahedron': 6,
    'Prism': 6,
}

GDICT2 = {
    'Point': 1,
    'Line': 4,
    'Triangle': 7,
    'Quadrangle': 7,
    'Tetrahedron': 15,
    'Hexahedron': 34,
    'Prism': 21,
}


MODEL_NAME = 'mspa'


def add_integration(integration, name, group_dict, itype='Gauss'):
    i0 = integration.add(name)
    item = i0.add()
    item_case = item.add(itype)
    ici = item_case.add()
    for element, value in group_dict.items():
        ici.add(GeoElement=element, NumberOfPoints=value)


def setup_planes():
    name = 'CutPlane'
    plugin.set_number(name, 'A', 0.0)
    plugin.set_number(name, 'B', 0.0)
    plugin.set_number(name, 'C', 1.0)
    plugin.set_number(name, 'D', -0.0145)
    plugin.set_number(name, 'View', 0)
    plugin.run(name)
    # name = 'ModulusPhase'
    # plugin.set_number(name, 'RealPart', 0)
    # plugin.set_number(name, 'ImaginaryPart', 1)
    # plugin.set_number(name, 'View', 2)
    # plugin.run(name)
    option.set_string('View[2].Name', 'e_amp')
    # option.set_number('View[2].ScaleType', 2)
    option.set_number('View[2].ForceNumComponents', 9)


def setup_onelab():
    gmsh.initialize()
    onelab.set(
        """
        [
            {
                "type": "number",
                "name": "Model/Frequency",
                "values": [137.1],
                "min": 141.5,
                "max": 141.0,
                "step": 10.0,
                "index": 0,
                "clients": {"Gmsh": 0}
            },
            {
                "type": "number",
                "name": "Model/epr",
                "values": [1.05],
                "min": 1.0,
                "max": 1.5,
                "step": 0.1,
                "index": 0,
                "clients": {"Gmsh": 0}
            },
            {
                "type": "number",
                "name": "Model/Lambda",
                "readOnly": true,
                "index": 1,
                "clients": {"Gmsh": 0}
            },
            {
                "type": "number",
                "name": "Model/WaveNumber",
                "label": "Wave Number",
                "readOnly": true,
                "index": 2,
                "clients": {"Gmsh": 0}
            }
        ]
        """
    )


def setup_plugins(r, wavenumber):
    name = 'CutSphere'
    plugin.set_number(name, 'Xc', 0.0)
    plugin.set_number(name, 'Yc', 0.0)
    plugin.set_number(name, 'Xc', 0.0)
    plugin.set_number(name, 'R', r)

    plugin.set_number(name, 'View', 0)
    plugin.run(name)
    plugin.set_number(name, 'View', 1)
    plugin.run(name)

    name = 'NearToFarField'
    plugin.set_number(name, 'Wavenumber', wavenumber)
    plugin.set_number(name, 'RFar', 1)
    plugin.set_number(name, 'NumPointsPhi', 60)  # 50
    plugin.set_number(name, 'NumPointsTheta', 30)  # 25
    plugin.set_number(name, 'EView', 3)
    plugin.set_number(name, 'HView', 4)
    plugin.set_number(name, 'Normalize', 0)
    plugin.set_number(name, 'dB', 1)
    plugin.run(name)


setup_onelab()
antenna = Mspa(MODEL_NAME)

pro = Problem()
pro.filename = MODEL_NAME + '.pro'
pro.include('defines.pro')

groups = model.get_physical_groups()
for g in groups:
    tag = g[1]
    name = model.get_physical_name(g[0], g[1])
    pro.group.add(name, tag)
pro.group.Region('SurBC', 'SkinFeed')
pro.group.Region(
    'DomainTot', ['Substrate', 'SkinFeed', 'Air', 'Pml', 'SigmaInf'])
pro.group.Region('DomainC', [])  # TODO remove
pro.group.Region('Domain', ['Substrate', 'Air', 'Pml'])
pro.group.define('DomainS')  # TODO remove
pro.group.define('SurS')  # TODO remove
pro.group.ElementsOf('TrGr', 'Domain', OnOneSideOf='SkinFeed')


fvar = {}
fvar['mu0'] = mu_0
fvar['nu0'] = 1.0 / mu_0
fvar['ep0'] = epsilon_0
# fvar['epr'] = 1.5  # 1.5  # Dielectric constant for FR4 is ~4.5

box = model.occ.get_bounding_box(*antenna.tags['vol_air'])

fvar['pml_xmax'] = box[3]
fvar['pml_ymax'] = box[4]
fvar['pml_zmax'] = box[5]
fvar['pml_xmin'] = box[0]
fvar['pml_ymin'] = box[1]
fvar['pml_zmin'] = box[2]
dc = 0.0  # 0.035e-3
gap = antenna.dims['gap']
fvar['gap'] = gap  # TODO refactor it
fvar['pml_delta'] = 0.2
fvar['air_boundary'] = 1.3
fvar['zl'] = 50.0  # Ohm load resistance

f = pro.function

for name, value in fvar.items():
    f.constant(name, value)

f.add('I', f.Complex(0.0, 1.0))
f.add('epsilon', 'ep0', region=['Air', 'SkinFeed', 'SigmaInf'])
f.add('epsilon', 'epr * ep0', region=['Substrate'])
f.add('nu', 'nu0', region=['Air', 'Substrate', 'SkinFeed', 'SigmaInf'])

f.add('sigma', '6.0e7')  # Copper 6.0e7
f.define('js0')  # TODO remove
f.define('ks0')  # TODO remove
f.define('nxh')  # TODO remove

f.add('r', f.Sqrt('X[]^2 + Y[]^2 + Z[]^2'))
f.add('dumping_r', '(r[] >= air_boundary) ? 1.0 / (pml_delta - (r[] - air_boundary)) : 0.0')
f.add('cx', f.Complex(1.0, '-dumping_r[] / k0'))
f.add('cy', f.Complex(1.0, '-dumping_r[] / k0'))
f.add('cz', f.Complex(1.0, '-dumping_r[] / k0'))

f.add('tens', f.TensorDiag('cy[] * cz[] / cx[]',
                           'cx[] * cz[] / cy[]',
                           'cx[] * cy[] / cz[]'))
f.add('epsilon', 'ep0 * tens[]', region='Pml')
f.add('nu', 'nu0 / tens[]', region='Pml')

y_feed = antenna.dims['d_feed']  # - 0.5 * antenna.dims['w_path']
f.constant('y_feed', y_feed)

# f.add('r_xy', f.Sqrt('X[]^2 + (Y[] + y_feed)^2'))
# f.add('BC_Fct_e', f.Vector('X[] / r_xy[] / gap',
#       '(Y[] + y_feed) / r_xy[] / gap', 0.0))

# f.add('dr', f.Vector('-(Y[] + y_feed) / r_xy[] / gap',
#       'X[] / r_xy[] / gap', 0.0), region=['SkinFeed'])

f.add('r_xy', f.Sqrt('(X[] + y_feed)^2 + (Y[] + y_feed)^2'))
f.add('BC_Fct_e', f.Vector('(X[] + y_feed) / r_xy[] / gap',
      '(Y[] + y_feed) / r_xy[] / gap', 0.0))

f.add('dr', f.Vector('-(Y[] + y_feed) / r_xy[] / gap',
      '(X[] + y_feed) / r_xy[] / gap', 0.0), region=['SkinFeed'])

constr = pro.constraint
ef = constr.add('ElectricField')
c0 = ef.add()
c0.add(Region='SkinFeed', Type='AssignFromResolution',
       NameOfResolution='Microwave_e_BC')
c0.add(Region='SkinConductor', Type='Assign', Value=0.0)
c0.add(Region='SigmaInf', Type='Assign', Value=0.0)

jacobian = pro.jacobian
for js, s in enumerate(['Vol', 'Sur']):
    jacobian.add(Name=('J' + s))
    jacobian.items[js].add()
    jacobian.items[js].cases[0].add(Region="All", Jacobian=s)

fspace = pro.functionspace
fs = fspace.add('Hcurl_e', Type='Form1')
fs.add_basis_function(
    Name='se',
    NameOfCoef='ee',
    Function='BF_Edge',
    Support='DomainTot',
    Entity='EdgesOf[All]'
)
fs.add_constraint(NameOfCoef='ee', EntityType='EdgesOf',
                  NameOfConstraint='ElectricField')
fs = fspace.add('Hcurl_h', Type='Form1')
fs.add_basis_function(
    Name='sh',
    NameOfCoef='he',
    Function='BF_Edge',
    Support='DomainTot',
    Entity='EdgesOf[All]'
)

add_integration(pro.integration, 'I1', GDICT1)
add_integration(pro.integration, 'I2', GDICT2)

formulation = pro.formulation
f = formulation.add('Microwave_e_BC', Type='FemEquation')
q = f.add_quantity()
q.add(Name='e', Type='Local', NameOfSpace='Hcurl_e')
e = f.add_equation()
e.add('Galerkin', '',
      'Dof{e} , {e}',
      In='SurBC', Integration='I2', Jacobian='JSur')
e.add('Galerkin', '',
      '-BC_Fct_e[] , {e}',
      In='SurBC', Integration='I2', Jacobian='JSur')

f = formulation.add('Microwave_e', Type='FemEquation')
q = f.add_quantity()
q.add(Name='e', Type='Local', NameOfSpace='Hcurl_e')
q.add(Name='h', Type='Local', NameOfSpace='Hcurl_h')

e = f.add_equation()
e.add('Galerkin', '',
      'nu[] * Dof{d e} , {d e}',
      In='Domain', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDof',
      'sigma[] * Dof{e}, {e}',
      In='DomainC', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDtDof',
      'epsilon[] * Dof{e} , {e}',
      In='Domain', Integration='I1', Jacobian='JVol')
e.add('Galerkin', '',
      'Dof{h} , {h}',
      In='TrGr', Integration='I1', Jacobian='JVol')
e.add('Galerkin', '',
      '-I[] * nu[] * Dof{d e} / (2.0 * Pi * freq), {h}',
      In='TrGr', Integration='I1', Jacobian='JVol')

resolution = pro.resolution

res = resolution.add('Microwave_e_BC; Hidden 1')
s = res.add_system()
s.add(Name='B', NameOfFormulation='Microwave_e_BC', DestinationSystem='A')

operation = res.add_operation()
operation.Generate('B')
operation.Solve('B')
operation.TransferSolution('B')

res = resolution.add('Analysis')
s = res.add_system()
s.add(Name='A', NameOfFormulation='Microwave_e',
      Type='Complex', Frequency='freq')

operation = res.add_operation()
operation.CreateDirectory('build')
operation.Generate('A')
operation.Solve('A')
operation.SaveSolution('A')

pp = pro.postprocessing
ppi = pp.add('Microwave_e', 'Microwave_e')
quantity = ppi.add()
quantity.add(Name='e', Type='Local',
             Value='{e}',
             In='DomainTot', Jacobian='JVol')
quantity.add(Name='h', Type='Local',
             Value='I[] * nu[] * {d e} / (2.0 * Pi * freq)',
             In='Domain', Jacobian='JVol')

# quantity.add(Name='e_norm', Type='Local',
#              Value='Norm[{e}]',
#              In='Domain', Jacobian='JVol')
# admittance
quantity.add(Name='y', Type='Integral',
             Value='{h} * dr[]', In='SkinFeed',
             Jacobian='JSur', Integration='I2')
quantity.add(Name='s11', Type='Term',
             Value='20.0 * Log10[Norm[(1.0 - zl * $y) / (1.0 + zl * $y)]]', In='SkinFeed')

po = pro.postoperation
poi = po.add('Microwave_e', 'Microwave_e')
poi0 = poi.add()
poi0.add('e', OnElementsOf='Region[{Domain}]', File='./build/e.pos')  # , -Pml
poi0.add('h', OnElementsOf='Region[{Domain}]', File='./build/h.pos')  # , -Pml
# poi0.add('e', OnLine='{{0.0, 0.0, 0.02} {0.0, 0.0, 1.1}} {100}',
#          File='./build/e_linez.pos')
# poi0.add('h', OnLine='{{0.0, 0.0, 0.02} {0.0, 0.0, 1.1}} {100}',
#          File='./build/h_linez.pos')
# poi0.add('e', OnLine='{{0.0, 0.0, 0.2} {0.0, 0.0, 1.2}} {100}', Format='SimpleTable',
#          File='./build/e_linez.txt')
# poi0.add(
#     'e', OnSection='{{0.0, 0.0, 0.0} {1.0, 0.0, 0.0} {0.0, 1.0, 0.0}}', File='./build/e_norm.pos')
poi0.add('y[SkinFeed]', OnGlobal='', Format='FrequencyTable',
         StoreInVariable='$y', File='./build/y.txt')
poi0.add('s11', OnRegion='SkinFeed', Format='FrequencyTable',
         StoreInVariable='$s11', SendToServer='"s11"', File='./build/s11.txt')

pro.make_file()
pro.write_file()
gmsh.open(pro.filename)
model.set_current(MODEL_NAME)

model.mesh.generate(3)
gmsh.write(f'{MODEL_NAME}.msh')
onelab.run()
setup_planes()
setup_plugins(1.1, onelab.get_number('Model/WaveNumber')[0])


# result = []
# for s in np.arange(90.0, 200.0, 10.0):
#     # onelab.set_number('Model/Frequency', [s])
#     # antenna.refresh()
#     antenna.d_feed = s * 0.001
#     model.mesh.generate(3)
#     gmsh.write(f'{MODEL_NAME}.msh')
#     model.set_current(MODEL_NAME)
#     onelab.run()
#     s11 = onelab.get_number('s11')[0]
#     size = onelab.get_number('Model/FeedDistance')[0]
#     freq = onelab.get_number('Model/Frequency')[0]
#     result.append([freq, size, s11])
#     # d = np.loadtxt('/home/taran/work/gmsh/mspa/build/e_linez.txt')
#     # np.savetxt(f'/home/taran/work/gmsh/mspa/build/e_linez_{int(s):03d}.txt', d)
#     # v1 = np.vectorize(complex)(d[:, 3], d[:, 4])
#     # v2 = np.vectorize(complex)(d[:, 6], d[:, 7])
#     # angle = np.degrees(np.angle(-v1 * v2))
#     # v1 = np.vectorize(complex)(d[:, 3], d[:, 4])
#     # v2 = np.vectorize(complex)(d[:, 6], d[:, 7])
#     # v = -v1 * v2
#     # angle = np.degrees(np.angle(v2) - np.angle(v1))
#     # angle[angle < 0.0] += 360.0
#     # pyplot.plot(angle, 'o-')
#     # pyplot.grid()
#     # pyplot.savefig(f'v1v2_{int(s):03d}.png')
#     # pyplot.cla()
#     # pyplot.clf()
#     # pyplot.plot(v, 'o-')
#     # pyplot.grid()
#     # pyplot.savefig(f'v_{int(s):03d}.png')
#     # pyplot.cla()
#     # pyplot.clf()
#     # result.append((size, angle))
# pprint(result)
# gmsh.finalize()
# exit(0)


# def check_event():
#     action = onelab.get_string('ONELAB/Action')
#     if len(action) == 0:
#         pass
#     elif action[0] == 'check':
#         # onelab.clear('ONELAB/Action')
#         # createGeometryAndMesh()
#         gmsh.graphics.draw()
#     elif action[0] == 'compute':
#         print(f'{action[0]}')
#         onelab.run()
#     onelab.clear('ONELAB/Action')
#     return True


if "-nopopup" not in sys.argv:
    gmsh.fltk.initialize()
    gmsh.fltk.run()
    # while gmsh.fltk.isAvailable() and check_event():
    #     gmsh.fltk.wait()

gmsh.finalize()
