import os
import sys
from pygetdp import Group, Function, Problem
from pygetdp.helpers import build_example_png, print_html
from math import pi
import gmsh
import mspa


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


def add_integration(integration, name, group_dict, itype='Gauss'):
    i0 = integration.add(name)
    item = i0.add()
    item_case = item.add(itype)
    ici = item_case.add()
    for element, value in group_dict.items():
        ici.add(GeoElement=element, NumberOfPoints=value)


pro = Problem()
pro.filename = mspa.MODEL_NAME + '.pro'
groups = gmsh.model.getPhysicalGroups()
for g in groups:
    tag = g[1]
    name = gmsh.model.getPhysicalName(g[0], g[1])
    pro.group.add(name, tag)
pro.group.Region('SurBC', 'SkinFeed')
pro.group.Region('DomainTot', ['SkinFeed', 'Air', 'Pml'])
pro.group.Region('DomainC', [])  # TODO remove
pro.group.Region('Domain', ['Substrate', 'Air', 'Pml'])
pro.group.define('DomainS')  # TODO remove
pro.group.define('SurS')  # TODO remove
pro.group.ElementsOf('TrGr', 'Domain', OnOneSideOf='SkinFeed')

freq = 1.6e9  # Hz~1.0/s
cvel = 299792458.0  # m/s

fvar = {}
fvar['mu0'] = 4.0e-7 * pi
fvar['nu0'] = 1.0 / fvar['mu0']
fvar['ep0'] = 8.854187817e-12
fvar['epr'] = 4.7  # Dielectric constant for FR4 is ~4.5
fvar['zl'] = 50.0  # Ohms load resistance
fvar['eta0'] = 120.0 * pi  # eta0 = Sqrt(mu0/eps0)
fvar['freq'] = freq
fvar['k0'] = 2.0 * pi * freq / cvel

fvar['pml_xmax'] = 0.055
fvar['pml_ymax'] = 0.042
fvar['pml_zmax'] = 0.010
fvar['pml_xmin'] = -0.003
fvar['pml_ymin'] = -0.003
fvar['pml_zmin'] = -0.003
fvar['pml_delta'] = 12.0 * (mspa.dc + mspa.dh)


f = pro.function

for name, value in fvar.items():
    f.constant(name, value)

f.add('I', f.Complex(0, 1))
f.add('epsilon', 'ep0', region=['Air', 'SkinFeed'])
f.add('epsilon', 'epr * ep0', region='Substrate')
f.add('nu', 'nu0', region=['Air', 'Substrate', 'SkinFeed'])

f.add('sigma', '6.0e7')  # Copper
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
                           'cx[] * cz[] / cy[]', 'cx[] * cy[] / cz[]'))
f.add('epsilon', 'ep0 * tens[]', region='Pml')
f.add('nu', 'nu0 / tens[]', region='Pml')

# D5 = 1.40 * mm ;
# V0 = 1 ; delta_gap = D5 ;
# BC_Fct_e[] =  V0/delta_gap * Vector[1, 0, 0] ;
f.add('BC_Fct_e', '1.0 / 1.4e-3 * ' + f.Vector(1.0, 0.0, 0.0))
# print_html(f.code)
# exit(0)

constr = pro.constraint
ef = constr.add('ElectricField')
c0 = ef.add()
c0.add(Region='SkinFeed', Type='AssignFromResolution',
       NameOfResolution='Microwave_e_BC')
c0.add(Region='Conductor', Type='Assign', Value=0.0)
c0.add(Region='Pml', Type='Assign', Value=0.0)

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
e.add('Galerkin', 'DtDof', 'nu[] * Dof{d e} , {d e}', In='Domain',
      Integration='I1', Jacobian='JVol')
e.add('Galerkin', '', 'sigma[] * Dof{e}, {e}',
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
e.add('Galerkin', '', '-I[]*nu[]*Dof{d e}/(2*Pi*freq), {h}',
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

# pp = pro.postprocessing
# pp.add('pp_test', 'Microwave_e')
# po = pro.postoperation
# poi = po.add('po_test', 'pp_test')
# poi.add().add("v", OnElementsOf="Domain", File="test.pos")


pro.make_file()
pro.write_file()
# gmsh.open(pro.filename)
# gmsh.onelab.run()
# gmsh.model.setCurrent('mspa.geo')


if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
