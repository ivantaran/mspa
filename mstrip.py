
from pygetdp import Group, Function, Problem
from pygetdp.helpers import build_example_png, print_html
from scipy.constants import mu_0, epsilon_0, pi, speed_of_light
import gmsh
import mspa
import numpy as np
import os
import sys

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


def _setup_plugins(box, wavenumber):
    p = gmsh.plugin

    name = 'CutBox'
    p.setNumber(name, 'NumPointsU', 100)
    p.setNumber(name, 'NumPointsV', 100)
    p.setNumber(name, 'NumPointsW', 5)

    xmin = box[3]
    ymin = box[4]
    zmin = box[5]
    xmax = box[0]
    ymax = box[1]
    zmax = box[2]

    p.setNumber(name, 'X0', xmin)
    p.setNumber(name, 'Y0', ymin)
    p.setNumber(name, 'Z0', zmin)

    p.setNumber(name, 'X1', xmax)
    p.setNumber(name, 'Y1', ymin)
    p.setNumber(name, 'Z1', zmin)

    p.setNumber(name, 'X2', xmin)
    p.setNumber(name, 'Y2', ymax)
    p.setNumber(name, 'Z2', zmin)

    p.setNumber(name, 'X3', xmin)
    p.setNumber(name, 'Y3', ymin)
    p.setNumber(name, 'Z3', zmax)

    p.setNumber(name, 'ConnectPoints', 1)
    p.setNumber(name, 'Boundary', 1)

    p.setNumber(name, 'View', 0)
    p.run(name)
    p.setNumber(name, 'View', 1)
    p.run(name)

    name = 'NearToFarField'
    p.setNumber(name, 'Wavenumber', wavenumber)
    p.setNumber(name, 'RFar', 1)
    p.setNumber(name, 'NumPointsPhi', 30)  # 50
    p.setNumber(name, 'NumPointsTheta', 15)  # 25
    p.setNumber(name, 'EView', 2)
    p.setNumber(name, 'HView', 3)
    p.setNumber(name, 'Normalize', 1)
    p.setNumber(name, 'dB', 1)
    p.run(name)


model = mspa.Mspa(MODEL_NAME)

pro = Problem()
pro.filename = MODEL_NAME + '.pro'
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

freq = 1.480e9  # 1.575e9

fvar = {}
fvar['mu0'] = mu_0
fvar['nu0'] = 1.0 / fvar['mu0']
fvar['ep0'] = epsilon_0
fvar['epr'] = 3.38  # Dielectric constant for FR4 is ~4.5
fvar['freq'] = freq
fvar['k0'] = 2.0 * pi * freq / speed_of_light

box = gmsh.model.occ.getBoundingBox(*model.tags['vol_air'])

fvar['pml_xmax'] = box[3]
fvar['pml_ymax'] = box[4]
fvar['pml_zmax'] = box[5]
fvar['pml_xmin'] = box[0]
fvar['pml_ymin'] = box[1]
fvar['pml_zmin'] = box[2]
dc = 0.0  # 0.035e-3
dh = model.dims['d']
fvar['dh'] = dh
fvar['pml_delta'] = 0.05  # 12.0 * (dc * 2.0 + dh)


f = pro.function

for name, value in fvar.items():
    f.constant(name, value)

f.add('I', f.Complex(0, 1))
f.add('epsilon', 'ep0', region=['Air', 'SkinFeed', 'SigmaInf'])
f.add('epsilon', 'epr * ep0', region='Substrate')
f.add('nu', 'nu0', region=['Air', 'Substrate', 'SkinFeed', 'SigmaInf'])

f.add('sigma', '6.0e7')  # Copper 6.0e7
f.define('js0')  # TODO remove
f.define('ks0')  # TODO remove
f.define('nxh')  # TODO remove

f.add('DampingProfileX', '(X[] >= pml_xmax || X[] <= pml_xmin) ? (X[] >= pml_xmax ? 1.0 / (pml_delta - (X[] - pml_xmax)): 1.0 / (pml_delta - (pml_xmin - X[]))): 0.0')
f.add('DampingProfileY', '(Y[] >= pml_ymax || Y[] <= pml_ymin) ? (Y[] >= pml_ymax ? 1.0 / (pml_delta - (Y[] - pml_ymax)): 1.0 / (pml_delta - (pml_ymin - Y[]))): 0.0')
f.add('DampingProfileZ', '(Z[] >= pml_zmax || Z[] <= pml_zmin) ? (Z[] >= pml_zmax ? 1.0 / (pml_delta - (Z[] - pml_zmax)): 1.0 / (pml_delta - (pml_zmin - Z[]))): 0.0')

f.add('cx', f.Complex(1.0, '-DampingProfileX[] / k0'))
f.add('cy', f.Complex(1.0, '-DampingProfileY[] / k0'))
f.add('cz', f.Complex(1.0, '-DampingProfileZ[] / k0'))
f.add('tens', f.TensorDiag('cy[] * cz[] / cx[]',
                           'cx[] * cz[] / cy[]',
                           'cx[] * cy[] / cz[]'))
f.add('epsilon', 'ep0 * tens[]', region='Pml')
f.add('nu', 'nu0 / tens[]', region='Pml')
f.add('BC_Fct_e', f.Vector(0.0, 0.0, 1.0 / dh))

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
e.add('Galerkin', '', 'Dof{e} , {e}', In='SurBC',
      Integration='I2', Jacobian='JSur')
e.add('Galerkin', '', '-BC_Fct_e[] , {e}',
      In='SurBC', Integration='I2', Jacobian='JSur')

f = formulation.add('Microwave_e', Type='FemEquation')
q = f.add_quantity()
q.add(Name='e', Type='Local', NameOfSpace='Hcurl_e')
q.add(Name='h', Type='Local', NameOfSpace='Hcurl_h')

e = f.add_equation()
e.add('Galerkin', '', 'nu[] * Dof{d e} , {d e}', In='Domain',
      Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDof', 'sigma[] * Dof{e}, {e}',
      In='DomainC', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDtDof', 'epsilon[] * Dof{e} , {e}',
      In='Domain', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDof', 'js0[], {e}',
      In='DomainS', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDof', '-ks0[] , {d e}',
      In='DomainS', Integration='I1', Jacobian='JVol')
e.add('Galerkin', 'DtDof', '-nxh[] , {e}',
      In='SurS', Integration='I1', Jacobian='JSur')

# // store magnetic field for Admitance computation (Yin)
e.add('Galerkin', '', 'Dof{h} , {h}', In='TrGr',
      Integration='I1', Jacobian='JVol')
e.add('Galerkin', '', '-I[] * nu[] * Dof{d e} / (2.0 * Pi * freq), {h}',
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
# operation.PostOperation('Microwave_e')

pp = pro.postprocessing
ppi = pp.add('Microwave_e', 'Microwave_e')
quantity = ppi.add()
quantity.add(Name='e', Type='Local',
             Value='{e}', In='DomainTot', Jacobian='JVol')
quantity.add(Name='h_from_e', Type='Local',
             Value='I[] * nu[] * {d e} / (2.0 * Pi * freq)', In='Domain', Jacobian='JVol')
# quantity.add(Name='exh', Type='Local',
#              Value='CrossProduct[{e}, Conj[I[] * nu[] * {d e} / (2.0 * Pi * freq)]]',
#              In='Domain', Jacobian='JVol')


po = pro.postoperation
poi = po.add('Microwave_e', 'Microwave_e')
poi0 = poi.add()
poi0.add('e', OnElementsOf='Region[{Domain, -Pml}]', File='./build/e.pos')
poi0.add(
    'h_from_e', OnElementsOf='Region[{Domain, -Pml}]', File='./build/h_pml.pos')
# poi0.add('exh', OnElementsOf='Region[{Domain, -Pml}]',
#          File='./build/exh_pml.pos')


gmsh.model.mesh.generate(3)
gmsh.write(MODEL_NAME + '.msh')
pro.make_file()
pro.write_file()
gmsh.open(pro.filename)

# gmsh.model.mesh.generate(1)
# gmsh.merge('./build/e.pos')
# gmsh.merge('./build/h_pml.pos')

gmsh.onelab.run()
gmsh.model.setCurrent(MODEL_NAME)

minimal_box = True
if minimal_box:
    box = gmsh.model.occ.getBoundingBox(
        *model.tags['vol_substrate'])  # sur_patch, vol_substrate
else:
    airbox = gmsh.model.occ.getBoundingBox(*model.tags['vol_substrate'])
    box = [0.0] * 6
    eps = 1.0e-5
    for i in range(3):
        box[i] = airbox[i] - eps
        box[i + 3] = airbox[i + 3] + eps
    # box[2] += 1.0e-3
    # box[5] += 1.0e-3
_setup_plugins(box, fvar['k0'])

# gmsh.onelab.run()
if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
