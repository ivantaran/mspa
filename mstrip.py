import os
import sys
from pygetdp import Group, Function, Problem
from pygetdp.helpers import build_example_png, print_html
from math import pi
import gmsh
import mspa

pro = Problem()
pro.filename = mspa.MODEL_NAME + '.pro'
groups = gmsh.model.getPhysicalGroups()
for g in groups:
    tag = g[1]
    name = gmsh.model.getPhysicalName(g[0], g[1])
    pro.group.add(name, tag)

# print_html(group.code)

fvar = {}
fvar['mu0'] = 4.0e-7 * pi
fvar['nu0'] = 1.0 / fvar['mu0']
fvar['ep0'] = 8.854187817e-12
fvar['epr'] = 4.7  # Dielectric constant for FR4 is ~4.5
fvar['zl'] = 50.0  # Ohms load resistance
fvar['eta0'] = 120.0 * pi  # eta0 = Sqrt(mu0/eps0)
fvar['freq'] = 1.6e9  # MHz

f = pro.function

for name, value in fvar.items():
    f.constant(name, value)

f.add('I', f.Complex(0, 1))
f.add('epsilon', 'ep0', region=['Air', 'SkinFeed', 'Pml'])
f.add('epsilon', 'epr * ep0', region='Substrate')

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
fs0 = fspace.add('Hcurl_e', Type='Form1')
fs0.add_basis_function(
    Name='se',
    NameOfCoef='ee',
    Function='BF_Edge',
    Support='DomainTot',
    Entity='EdgesOf[All]'
)
fs0.add_constraint(NameOfCoef='ee', EntityType='EdgesOf',
                   NameOfConstraint='ElectricField')

integration = pro.integration
i0 = integration.add("I2")
item = i0.add()
item_case = item.add("Gauss")
ici = item_case.add()

gdict = {
    'Point': 1,
    'Line': 4,
    'Triangle': 7,
    'Quadrangle': 7,
    'Tetrahedron': 15,
    'Hexahedron': 34,
    'Prism': 21,
}

for name, value in gdict.items():
    ici.add(GeoElement=name, NumberOfPoints=value)

formulation = pro.formulation
f0 = formulation.add('Microwave_e_BC', Type='FemEquation')
q = f0.add_quantity()
q.add(Name='e', Type='Local', NameOfSpace='Hcurl_e')
e = f0.add_equation()
e.add('Galerkin', 'Dof{e} , {e}', In='SkinFeed',
      Integration='I2', Jacobian='JSur')
e.add('Galerkin', '-BC_Fct_e[] , {e}',
      In='SkinFeed', Integration='I2', Jacobian='JSur')

resolution = pro.resolution
res0 = resolution.add('Microwave_e_BC')
s = res0.add_system()
# s.add(Name='B', NameOfFormulation='Microwave_e_BC', DestinationSystem='A')
s.add(Name='B', NameOfFormulation='Microwave_e_BC')
operation = res0.add_operation()
operation.Generate('B')
operation.Solve('B')
operation.SaveSolution('B')
# operation.TransferSolution('B')

pro.make_file()
pro.write_file()
# gmsh.merge(pro.filename)

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
