from ..TB.latticeTB import *
from plotTB import *
from eigTB import *
from math import sqrt


# dimer chain
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [1, 0]}]
prim_vec = {'norm': 2, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
n1, n2 = 6, 6
lat.get_lattice(n1=n1, n2=n2)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)


lat.get_lattice(n1=n1)
fig_lat = plot.lattice(ms=15, figsize=(10, 2))
save = saveTB(sys=sys, dir_name='chain')

# square lattice with 1 site per unit cell
unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
prim_vec = {'norm': 1, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='square1')
n1, n2 = 8, 8
lat.get_lattice(n1=n1, n2=n2)
fig_lat = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')

# face centered hexagon
unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
prim_vec = {'norm': 1, 'angle': 120}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='hexa_line_center')
n1, n2 = 7, 7
lat.get_lattice(n1=n1, n2=n2)
fig_lat = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')

# square lattice with 4 sites per unit cell
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [1, 0]},
                 {'tag': b'c', 'r0': [0, 1]},
                 {'tag': b'd', 'r0': [1, 1]}]
prim_vec = {'norm': 2, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='square4')

n1, n2 = 4, 4
lat.get_lattice(n1=n1, n2=n2)
fig_lat = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')

# Face-centered square
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [1, 1]}]
prim_vec = {'norm': 2, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='face_center_square')

n1, n2 = 6, 6
lat.get_lattice(n1=n1, n2=n2)
lat.remove_sites_x(xlim=[0, 10])
lat.remove_sites_y(ylim=[0, 10])
fig_lat = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')

# graphene
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='graphene')

n1, n2 = 6, 6
lat.get_lattice(n1=n1, n2=n2)
fig_lat_ind = plot.lattice(ms=20, plt_index=True)
lat.remove_sites([27, 32, 33, 38, 39, 44])
lat.remove_dangling()
fig_lat_vac = plot.lattice(ms=20)
save.fig_lat(fig_lat_ind, 'lattice_ind')
save.fig_lat(fig_lat_vac, 'lattice_index_vac')

# Lieb lattice
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [1, 0]},
                 {'tag': b'c', 'r0': [0, 1]}]
prim_vec = {'norm': 2, 'angle': 90}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='lieb')

n1, n2 = 5, 5
lat.get_lattice(n1=n1, n2=n2)
fig_lat = plot.lattice(ms=20)
lat.remove_dangling()
fig_lat_dang = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')
save.fig_lat(fig_lat_dang, 'lattice_dang')

# Kagome lattice
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [1, 0]},
                 {'tag': b'c', 'r0': [1/2., sqrt(3)/2]}]
prim_vec = {'norm': 2, 'angle': 60}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys)
save = saveTB(sys=sys, dir_name='kagome')
n1, n2 = 6, 5
lat.get_lattice(n1=n1, n2=n2)
fig_lat = plot.lattice(ms=20)
save.fig_lat(fig_lat, 'lattice')

# Line centered hexagon lattice
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                  {'tag': b'a', 'r0': [0.5*sqrt(3), 0.5]}, 
                  {'tag': b'c', 'r0': [0.25*sqrt(3), 0.25]},
                  {'tag': b'c', 'r0': [0.75*sqrt(3), 0.25]},
                  {'tag': b'c', 'r0': [0.5*sqrt(3), 1.]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
sys = eigTB(lat=lat)
plot = plotTB(sys=sys, colors=['b', 'b', 'r', 'r', 'r'])
save = saveTB(sys=sys, dir_name='graphene')

n1, n2 = 6, 6
lat.get_lattice(n1=n1, n2=n2)
lat.remove_dangling()
fig_lat_ind = plot.lattice(ms=15)
save.fig_lat(fig_lat_ind, 'lattice')
plt.show()