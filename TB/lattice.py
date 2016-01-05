import numpy as np
import matplotlib.pyplot as plt
import error_handling

PI = np.pi


class lattice():
    def __init__(self, unit_cell, prim_vec):
        '''
        Build up 1D or 2D lattice.
        Lattice is expressed by: 
        :math:`\mathbf{R} = n_1\mathbf{a}_1 + n_2\mathbf{a}_2`
        where :math:`\mathbf{a}_1` and `\mathbf{a}_2` are the two primitive vectors.
        where :math:`n_1` and `n_2` are the number of unit cells along
          `\mathbf{a}_1` and `\mathbf{a}_2`.
        :param unit_cell: List of dictionaries. One dictionary per site within the unit cell.
          Each dictionary has two keys: 
            * 'tag', label of the site.
            * 'r0', position. 
        :param prim_vec: Dictionary. Define the two primitive vectors.
           Dictionary with two keys: 
            * 'norm', norm of both primitive vectors: `|\mathbf{a}_1|=|\mathbf{a}_2| = norm` 
            * 'angle', angle between the two primitive:
                * :math:`\mathbf{a}_1| = norm (1, 0)` .
                * :math:`|\mathbf{a}_1|= norm (\cos angle, \sin \angle )`.
        '''
        error_handling.unit_cell(unit_cell)
        error_handling.prim_vec(prim_vec)
        self.unit_cell = unit_cell
        self.prim_vec = prim_vec
        self.tags = np.unique(np.array([dic['tag'] for dic in self.unit_cell]))
        self.n1, self.n2 = 0, 0  # number of unit cells along a1 and a2
        self.coor = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        self.sites = 0  # Site number

    def get_lattice(self, n1, n2=1):
        '''
        Get the lattice positions

        :param n1: Positive Integer. 
            Number of unit cells along :math:`\mathbf{a}_1`.
        :param n2: Positive Integer. Default value 1. 
            Number of unit cells along :math:`\mathbf{a}_2`.
        '''
        self.coor = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        error_handling.get_lattice(n1, n2)
        self.n1, self.n2 = n1, n2
        x = self.prim_vec['norm'] * np.arange(n1, dtype='f16')
        y = np.sin(PI / 180. * self.prim_vec['angle']) \
            * self.prim_vec['norm'] * np.arange(n2, dtype='f16')
        xx, yy = np.meshgrid(x, y)
        if self.prim_vec['angle'] == 0:
            xx += yy
        else:
            xx += yy / np.tan(PI / 180. * self.prim_vec['angle'])
        xx = np.ravel(xx)
        yy = np.ravel(yy)
        coor_tag = np.zeros(len(xx), dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        for dic in self.unit_cell:
            coor_tag['x'] = xx + dic['r0'][0]
            coor_tag['y'] = yy + dic['r0'][1]
            coor_tag['tag'] = dic['tag']
            self.coor = np.concatenate([self.coor, coor_tag])
        self.coor = np.sort(self.coor, order=('y', 'x'))        
        self.sites = len(self.coor['tag'])

    def add_sites(self, coor):
        '''
        Add sites.

        :param coor: Structured array with keys: {'x', 'y', 'tag'}.
        '''
        error_handling.coor(coor)
        self.coor = np.concatenate([self.coor, coor])
        self.sites += len(self.coor)

    def remove_sites(self, index):
        '''
        Remove sites defined by their indices
        (use method lattice(plt_index=True) of plotTB 
        to get access to the site indices).

        :param index: List of site indices to be removed.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.remove_sites(index, self.sites)
        mask = np.ones(self.sites, bool)
        mask[index] = False
        self.coor = self.coor[mask]
        self.sites = self.coor.size

    def remove_dangling(self):
        '''
        Remove dangling sites
        (sites connected with just another site).
        '''
        error_handling.empty_coor(self.coor)
        while True:
            dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
            dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
            dis = np.sqrt(dif_x ** 2 + dif_y ** 2)
            dis_unique = np.unique(dis)
            len_hop = dis_unique[1]
            ind = np.argwhere(np.isclose(dis, len_hop))
            dang = []
            for i in range(self.sites):
                if (ind[:, 0] == i).sum() == 1:
                    dang.append(i)
            self.coor = np.delete(self.coor, dang, axis=0)
            self.sites -= len(dang)
            if dang == []:
                break

    def shift_x(self, shift):
        '''
        Shift the x coordinates.

        :param shift: Real number. Shift value.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.real_number(shift, 'shift')
        self.coor['x'] += shift

    def shift_y(self, shift):
        '''
        Shift by *delta_x* the x coordinates.

        :param shift: Real number. Shift value.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.real_number(shift, 'shift')
        self.coor['y'] += shift

    def change_sign_x(self):
        '''
        Change x coordinates sign.
        '''
        error_handling.empty_coor(self.coor)
        self.coor['x'] *= -1

    def change_sign_y(self):
        '''
        Change y coordinates sign.
        '''
        error_handling.empty_coor(self.coor)
        self.coor['y'] *= -1

    def boundary_line(self, cx, cy, co):
        '''
        Select sites according to :math:`c_yy+c_xx > c_0`.

        :param cx: Real number. cx value.
        :param cy: Real number. cy value.
        :param co: Real number. co value.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.real_number(cx, 'cx')
        error_handling.real_number(cy, 'cy')
        error_handling.real_number(co, 'co')
        self.coor = self.coor[cy * self.coor['y'] + cx * self.coor['x'] > co]
        self.sites = len(self.coor)

    def ellipse_in(self, a, b):
        '''
        Select sites according to :math:`x^2/a^2+y^2/b^2 < 1`.

        :param a: Real number. a value.
        :param b: Real number. b value.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.positive_real(a, 'a')
        error_handling.positive_real(b, 'b')
        self.coor = self.coor[(self.coor['x'] / a) ** 2 +  (self.coor['y'] / b) ** 2 < 1.]
        self.sites = len(self.coor)

    def ellipse_out(self, a, b):
        '''
        Select sites according to :math:`x^2/a^2+y^2/b^2 > 1`.

        :param a: Real number. a value.
        :param b: Real number. b value.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.positive_real(a, 'a')
        error_handling.positive_real(b, 'b')
        self.coor = self.coor[(self.coor['x'] / a) ** 2 +  (self.coor['y'] / b) ** 2 >= 1.]
        self.sites = len(self.coor)

    def center(self):
        '''
        Fix the center of mass of the lattice at (0, 0).
        '''
        error_handling.empty_coor(self.coor)
        self.coor['x'] -= np.mean(self.coor['x'])
        self.coor['y'] -= np.mean(self.coor['y'])

    def __add__(self, other):
        '''
        Overloading operator +.
        '''
        error_handling.lat(other)
        error_handling.empty_coor(self.coor)
        error_handling.empty_coor(other.coor)
        sites = self.sites + other.sites
        tags = np.unique([self.tags, other.tags])
        coor = np.concatenate([self.coor, other.coor])
        lat = lattice(unit_cell=self.unit_cell, prim_vec=self.prim_vec)
        lat.add_sites(coor)
        lat.sites = self.sites
        lat.tags = self.tags
        return lat

    def __iadd__(self, other):
        '''
        Overloading operator +=.
        '''
        error_handling.lat(other)
        error_handling.empty_coor(self.coor)
        error_handling.empty_coor(other.coor)
        self.coor = np.concatenate([self.coor, other.coor])
        self.sites += other.sites
        self.tags = np.unique([self.tags, other.tags])
        return self

    def __sub__(self, other):
        '''
        Overloading operator -.
        '''
        error_handling.lat(other)
        error_handling.empty_coor(self.coor)
        error_handling.empty_coor(other.coor)
        tags = np.unique([self.tags, other.tags])
        boo = np.zeros(self.sites, bool)
        for i, c in enumerate(other.coor):
            boo += np.isclose(c['x'], self.coor['x']) & np.isclose(c['y'], self.coor['y'])
        coor = self.coor[np.logical_not(boo)]
        lat = lattice(unit_cell=self.unit_cell, prim_vec=self.prim_vec)
        lat.add_sites(coor)
        lat.sites = len(lat.coor)
        lat.tags = self.tags
        return lat

    def __isub__(self, other):
        '''
        Overloading operator -=.
        '''
        error_handling.lat(other)
        error_handling.empty_coor(self.coor)
        error_handling.empty_coor(other.coor)
        ind_remove = []
        boo = np.zeros(self.sites, bool)
        for i, c in enumerate(other.coor):
            boo += np.isclose(c['x'], self.coor['x']) & np.isclose(c['y'], self.coor['y'])
        self.coor = self.coor[np.logical_not(boo)]
        self.sites = sum(np.logical_not(boo))
        return self


"""

    def plot(self, ms=20, lw=5, fs=20, plt_index=None, figsize=None):
        '''
        Plot lattice
        '''
        # plot sites
        colors = ['b', 'r', 'g', 'y', 'm', 'k']
        fig, ax = plt.subplots(figsize=figsize)
        for color, dic in zip(colors, self.unit_cell):
            plt.plot(self.coor['x'][self.coor['tag'] == dic['tag']],
                        self.coor['y'][self.coor['tag'] == dic['tag']],
                       'o', color=color, ms=ms, markeredgecolor='none')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(self.coor['x'])-0.5, np.max(self.coor['x'])+0.5])
        ax.set_ylim([np.min(self.coor['y'])-0.5, np.max(self.coor['y'])+0.5])
        # plot indices
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sys.lat.sites)]
            for l, x, y in zip(indices, coor['x'], coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        #xa = ax.get_xaxis()
        #xa.set_major_locator(plt.MaxNLocator(integer=True))
        #ya = ax.get_yaxis()
        #ya.set_major_locator(plt.MaxNLocator(integer=True))
        plt.draw()
        return fig

from math import sin, sqrt, pi

'''
from math import sqrt

from eigTB import *
from plotTB import *


# ADD
n1, n2 = 5, 4
unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
prim_vec = {'norm': 1, 'angle': 0}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)

lat1.get_lattice(n1=n1, n2=1)
print(lat1.coor)
unit_cell = [{'tag': b'a', 'r0': [2, -4]}]
prim_vec = {'norm': 1, 'angle': 90}
lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
lat2.get_lattice(n1=1, n2=n2)
print(lat2.coor)
lat1 += lat2
lat1.plot()

'''
'''
# REMOVE
n1, n2 = 5, 2
unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
prim_vec = {'norm': 1, 'angle': 0}
lat3 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)

lat3.get_lattice(n1=n1, n2=1)
unit_cell = [{'tag': b'a', 'r0': [2, 0]}]
prim_vec = {'norm': 1, 'angle': 0}
lat4 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
lat4.get_lattice(n1=1, n2=n2)
lat3 -= lat4
print(lat3.coor)
lat3.plot()
'''
'''
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
n1, n2 = 7, 7
lat1.get_lattice(n1=n1, n2=n2)
fig_lat = lat1.plot(ms=20)


unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
n1, n2 = 3, 3
lat2.get_lattice(n1=n1, n2=n2)

#lat2.coor['x'] += sqrt(3)
#lat2.coor['x'] = np.round(lat2.coor['x'], 3)
lat2.coor['y'] += 3
fig_lat = lat2.plot(ms=20)

lat1 -= lat2
lat1.plot()


unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
'''
'''
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
n1, n2 = 2, 2
lat1.get_lattice(n1=n1, n2=n2)
#lat1.remove_dangling()
fig_lat = lat1.plot(ms=20)

unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                  {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
coor = np.array([(0, 0, 'a'), (1, 1, 'a'), (2, 2, 'a')], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])

lat1.add_sites(coor)
print(lat1.coor)
#lat1.remove_dangling()
fig_lat = lat1.plot(ms=20)
'''


from system import *
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)

n1, n2 = 10, 10
#sys1 = system(lat=lat1)

#sys1.set_hopping([{'n': 1, 't': 1}])
#lat1.remove_dangling()
#n = 21

# triangle zz
#lat1.get_lattice(n1=n+2, n2=n+2)
#lat1.boundary_line(cx=-sqrt(3), cy=-1, co=-3*(n+1))
#fig_lat = lat1.plot(ms=20)

# hexa zz
#lat1.get_lattice(n1=2*n, n2=2*n)
#lat1.boundary_line(cx=sqrt(3), cy=1, co=3*(n-1))
#lat1.boundary_line(cx=-sqrt(3), cy=-1, co=-9*n+1)

#tri ac


# hexa ac

#n2 = int(1 + (2*n * sqrt(3) - 2) / 2.5)
#lat1.get_lattice(n1=n, n2=10)



#square
#n2 = int(1 + (2*n * sqrt(3) - 2) / 2.5)
#lat1.get_lattice(n1=n, n2=10)
#(self.coor['x'] > DX * (n-1)-0.1) & (self.coor['x'] < DX * 2*n-0.1)
#lat1.boundary_line(cx=-1, cy=0, co=-2*sqrt(3)*n)
#lat1.boundary_line(cx=1, cy=0, co=sqrt(3)/2*(2*n-1)+.1)

#lat1.boundary_line(cx=-sqrt(3), cy=0, co=-9*n+1)
#self.sites = len(self.coor)
#ind = np.argwhere((self.coor['x'] > DX * (n-1)-0.1) & (self.coor['x'] < DX * 2*n-0.1))
#lat1.boundary_line(cx=-0.5, cy=+1, co=1)
#lat1.boundary_line(cx=1, cy=0, co=-10)
#lat1.boundary_line(cx=-1, cy=0, co=10)
#lat1.boundary_line(cx=0, cy=1, co=-5)
#lat1.boundary_line(cx=0, cy=-1, co=5)

#lat1.remove_dangling()

#circle
n = 3
#lat1.get_lattice(n1=2*n, n2=2*n)
#(self.coor['x'] > DX * (n-1)-0.1) & (self.coor['x'] < DX * 2*n-0.1)
#lat1.coor_center()
#lat1.ellipse_in(a=n*sqrt(3)/2, b=n*sqrt(3)/2)
#lat1.boundary_line(cx=1, cy=0, co=sqrt(3)/2*(2*n-1)+.1)
nn = 3 * n - 2
lat1.get_lattice(n1=2*nn, n2=2*nn)
'''
ind = np.argwhere((self.coor['x'] > DX * (2*nn-1)-.1) & 
                            (self.coor['y'] >  0.5 / DX  * self.coor['x']  - nn - .1) &
                            (self.coor['y'] < -0.5 / DX  * self.coor['x'] + 4 * nn - 1 + .1))
'''
lat1.boundary_line(cx=1, cy=0, co=sqrt(3)/2* (2*nn-1)-.1)
lat1.boundary_line(cx=-1/sqrt(3), cy=1, co=-nn -.1)
lat1.boundary_line(cx=-1/sqrt(3), cy=-1, co=-4*nn +1-.1)

fig_lat = lat1.plot(ms=20)
#lat1.line_lower(a=-0.5, b=4)
#fig_lat = lat1.plot(ms=20)

'''
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
n1, n2 = 30, 10
lat1.get_lattice(n1=n1, n2=n2)
lat1.coor_center()
fig_lat = lat1.plot(ms=20)
lat1.ellipse_in(a=12*sqrt(3)-0.5, b=11+0.5)
fig_lat = lat1.plot(ms=10)

lat1.remove_dangling()
'''
from math import sqrt
from plot import *
unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                 {'tag': b'b', 'r0': [0.5*sqrt(3), 0.5]}]
prim_vec = {'norm': sqrt(3), 'angle': 60}
lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
plt = plot(lat=lat)

n = 1 # no of plackets along one edge.
lat.get_lattice(n1=n+2, n2=n+2)

'''
"""
'''
class a:
    def __init__(self):
        pass

aa = a()
bb = a()
print(isinstance(aa, a))
print(isinstance(bb, a))
'''