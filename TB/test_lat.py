from lattice import *
from plot import *

n1, n2 = 5, 4
unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
prim_vec = {'norm': 1, 'angle': 0}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
lat1.get_lattice(n1=11, n2=1)
plt1 = plot(lat=lat1)
plt1.lattice()
plt1.show()