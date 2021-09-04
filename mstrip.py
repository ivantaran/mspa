
from numpy.lib.type_check import imag
from pygetdp import Group, Function, Problem
from pygetdp.helpers import build_example_png, print_html
from scipy.constants import mu_0, epsilon_0, pi, speed_of_light
import gmsh
# from comsol_patch_1575 import Mspa
from patch_137 import Mspa
import numpy as np
import os
import sys
from pprint import pprint
from matplotlib import pyplot

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


def _setup_planes():
    p = gmsh.plugin
    name = 'CutPlane'
    p.setNumber(name, 'A', 0.0)
    p.setNumber(name, 'B', 0.0)
    p.setNumber(name, 'C', 1.0)
    p.setNumber(name, 'D', 0.0)
    p.setNumber(name, 'View', 0)
    p.run(name)
    name = 'ModulusPhase'
    p.setNumber(name, 'RealPart', 0)
    p.setNumber(name, 'ImaginaryPart', 1)
    p.setNumber(name, 'View', 2)
    p.run(name)
    gmsh.option.setString('View[2].Name', 'e_amp')
    # gmsh.option.setNumber('View[2].ScaleType', 2)
    gmsh.option.setNumber('View[2].ForceNumComponents', 9)


def setup_onelab():
    gmsh.initialize()
    gmsh.onelab.set(
        """
        [
            {
                "type": "number",
                "name": "Model/Frequency",
                "values": [139.0],
                "min": 130.0,
                "max": 141.0,
                "step": 10.0,
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
            },
            {
                "type": "number",
                "name": "Model/CutRadius",
                "label": "Cut Radius",
                "values": [1.0],
                "min": 10.0,
                "max": 170.0,
                "step": 10.0,
                "clients": {"GetDP": 1}
            },
            {
                "type": "number",
                "name": "Model/FeedDistance",
                "label": "Feed Distance",
                "values": [120.0],
                "min": 50.0,
                "max": 400.0,
                "step": 10.0,
                "clients": {"GetDP": 1}
            },
            {
                "type": "number",
                "name": "Model/PatchSize",
                "label": "PatchSize",
                "values": [670.0],
                "min": 600.0,
                "max": 825.0,
                "step": 10.0,
                "index": 4,
                "clients": {"Gmsh": 0}
            }
        ]
        """
    )


def _setup_plugins(r, wavenumber):
    p = gmsh.plugin

    name = 'CutSphere'
    p.set_number(name, 'Xc', 0.0)
    p.set_number(name, 'Yc', 0.0)
    p.set_number(name, 'Xc', 0.0)
    p.set_number(name, 'R', r)

    p.set_number(name, 'View', 0)
    p.run(name)
    p.set_number(name, 'View', 1)
    p.run(name)

    name = 'NearToFarField'
    p.setNumber(name, 'Wavenumber', wavenumber)
    p.setNumber(name, 'RFar', 1)
    p.setNumber(name, 'NumPointsPhi', 30)  # 50
    p.setNumber(name, 'NumPointsTheta', 15)  # 25
    p.setNumber(name, 'EView', 2)
    p.setNumber(name, 'HView', 3)
    p.setNumber(name, 'Normalize', 0)
    p.setNumber(name, 'dB', 1)
    p.run(name)


setup_onelab()
model = Mspa(MODEL_NAME)

pro = Problem()
pro.filename = MODEL_NAME + '.pro'
pro.include('defines.pro')

groups = gmsh.model.getPhysicalGroups()
for g in groups:
    tag = g[1]
    name = gmsh.model.getPhysicalName(g[0], g[1])
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
fvar['epr'] = 1.5  # 1.5  # Dielectric constant for FR4 is ~4.5
# fvar['freq'] = freq
# fvar['k0'] = 2.0 * pi * freq / speed_of_light

box = gmsh.model.occ.getBoundingBox(*model.tags['vol_air'])

fvar['pml_xmax'] = box[3]
fvar['pml_ymax'] = box[4]
fvar['pml_zmax'] = box[5]
fvar['pml_xmin'] = box[0]
fvar['pml_ymin'] = box[1]
fvar['pml_zmin'] = box[2]
dc = 0.0  # 0.035e-3
gap = model.dims['gap']
fvar['gap'] = gap
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

y_feed = model.dims['d_feed'] - 0.5 * model.dims['w_path']
f.constant('y_feed', y_feed)

f.add('r_xy', f.Sqrt('X[]^2 + (Y[] + y_feed)^2'))
f.add('BC_Fct_e', f.Vector('X[] / r_xy[] / gap',
      '(Y[] + y_feed) / r_xy[] / gap', 0.0))
# f.add('BC_Fct_e', f.Vector(0.0, 0.0, 1.0 / gap))
f.add('dr', f.Vector('-(Y[] + y_feed) / r_xy[] / gap',
      'X[] / r_xy[] / gap', 0.0), region=['SkinFeed'])

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
# poi0.add('e', OnLine='{{0.0, 0.0, 0.2} {0.0, 0.0, 1.2}} {100}',
#          File='./build/e_linez.pos')
# poi0.add('e', OnLine='{{0.0, 0.0, 0.2} {0.0, 0.0, 1.2}} {100}', Format='SimpleTable',
#          File='./build/e_linez.txt')
poi0.add('y[SkinFeed]', OnGlobal='', Format='FrequencyTable',
         StoreInVariable='$y', File='./build/y.txt')
poi0.add('s11', OnRegion='SkinFeed', Format='FrequencyTable',
         StoreInVariable='$s11', SendToServer='"s11"', File='./build/s11.txt')

pro.make_file()
pro.write_file()
gmsh.open(pro.filename)
gmsh.model.set_current(MODEL_NAME)

# gmsh.model.mesh.generate(3)
# gmsh.write(f'{MODEL_NAME}.msh')
# gmsh.onelab.run()

# d = np.loadtxt('/home/taran/work/gmsh/mspa/build/e_linez.txt')
# v1 = np.vectorize(complex)(d[:, 3], d[:, 4])
# v2 = np.vectorize(complex)(d[:, 6], d[:, 7])
# v = -v1 * v2
# angle = np.degrees(np.angle(v2) - np.angle(v1))
# angle[angle < 0.0] += 360.0
# pyplot.plot(angle, 'o-')
# pyplot.plot(np.real(v1), np.imag(v1), 'o-')
# pyplot.plot(np.real(v2), np.imag(v2), 'o-')
# pyplot.xlim([-5, 5])
# pyplot.ylim([-5, 5])
# pyplot.grid()
# pyplot.show()

result = []
for s in np.arange(330.0, 150.0, -10.0):
    # gmsh.onelab.set_number('Model/CutRadius', [s])
    gmsh.onelab.set_number('Model/FeedDistance', [s])
    # gmsh.onelab.set_number('Model/Frequency', [s])
    # gmsh.onelab.set_number('Model/PatchSize', [s])
    model.refresh()
    gmsh.model.mesh.generate(3)
    gmsh.write(f'{MODEL_NAME}.msh')
    gmsh.model.set_current(MODEL_NAME)
    gmsh.onelab.run()
    s11 = gmsh.onelab.get_number('s11')[0]
    size = gmsh.onelab.get_number('Model/FeedDistance')[0]
    freq = gmsh.onelab.get_number('Model/Frequency')[0]
    result.append([freq, size, s11])
    # d = np.loadtxt('/home/taran/work/gmsh/mspa/build/e_linez.txt')
    # np.savetxt(f'/home/taran/work/gmsh/mspa/build/e_linez_{int(s):03d}.txt', d)
    # v1 = np.vectorize(complex)(d[:, 3], d[:, 4])
    # v2 = np.vectorize(complex)(d[:, 6], d[:, 7])
    # angle = np.degrees(np.angle(-v1 * v2))
    # v1 = np.vectorize(complex)(d[:, 3], d[:, 4])
    # v2 = np.vectorize(complex)(d[:, 6], d[:, 7])
    # v = -v1 * v2
    # angle = np.degrees(np.angle(v2) - np.angle(v1))
    # angle[angle < 0.0] += 360.0
    # pyplot.plot(angle, 'o-')
    # pyplot.grid()
    # pyplot.savefig(f'v1v2_{int(s):03d}.png')
    # pyplot.cla()
    # pyplot.clf()
    # pyplot.plot(v, 'o-')
    # pyplot.grid()
    # pyplot.savefig(f'v_{int(s):03d}.png')
    # pyplot.cla()
    # pyplot.clf()
    # result.append((size, angle))
pprint(result)
gmsh.finalize()
exit(0)

# _setup_planes()

# minimal_box = True
# if minimal_box:
#     box = gmsh.model.occ.getBoundingBox(
#         *model.tags['vol_substrate'])  # sur_patch, vol_substrate, vol_patch
# else:
#     airbox = gmsh.model.occ.getBoundingBox(*model.tags['vol_patch'])
#     box = [0.0] * 6
#     eps = 1.0e-3
#     for i in range(6):
#         box[i] = airbox[i]
#     box[2] = -eps * 10.0
#     box[5] = +eps * 10.0
_setup_plugins(1.1, gmsh.onelab.get_number('Model/WaveNumber')[0])


# def check_event():
#     action = gmsh.onelab.get_string('ONELAB/Action')
#     if len(action) == 0:
#         pass
#     elif action[0] == 'check':
#         # gmsh.onelab.clear('ONELAB/Action')
#         # createGeometryAndMesh()
#         gmsh.graphics.draw()
#     elif action[0] == 'compute':
#         print(f'{action[0]}')
#         gmsh.onelab.run()
#     gmsh.onelab.clear('ONELAB/Action')
#     return True


if "-nopopup" not in sys.argv:
    gmsh.fltk.initialize()
    gmsh.fltk.run()
    # while gmsh.fltk.isAvailable() and check_event():
    #     gmsh.fltk.wait()

gmsh.finalize()
