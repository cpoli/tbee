from lattice import *
from plot import *
from system import *
import numpy as np
from math import sqrt, pi

PI = np.pi
DX = 0.5 * np.sqrt(3)
DY = 0.5


def check_n(n):
    '''
    Test if the number *n* of plackets is a positive integer.

    :raises TypeError: Parameter *n* must be an integer.
    :raises ValueError: Parameter *n* must be a positive integer.
    '''
    if not isinstance(n, int):
        raise TypeError('\n\nParamer *n* must be an integer.')
    if n < 1:
        raise ValueError('\n\nParamer *n* must be a positive integer.')


def test_hop_linear_strain(t, beta):
    '''
    Test method *hop_linear_strain*.

    :raises TypeError: Parameter *t* must be an integer.
    :raises TypeError: Parameter *beta* must be a real number.
    '''
    if not isinstance(t, (int, float, complex)):
        raise TypeError('\n\nParamer *n* must be a number.')
    if not isinstance(beta, (int, float)):
        raise TypeError('\n\nParamer *n* must be a real number.')

def test_get_butterfly(t, N):
    '''
    Test method *hop_linear_strain*.

    :raises TypeError: Parameter *t* must be a real number.
    :raises TypeError: Parameter *N* must be an integer.
    :raises ValueError: Parameter *N* must be a positive integer.
    '''
    if not isinstance(t, (int, float)):
        raise TypeError('\n\nParamer *t* must be a number.')
    if not isinstance(N, int):
        raise TypeError('\n\nParamer *N* must be a real number.')
    if N < 1:
        raise ValueError('\n\nParamer *N* must be a positive integer.')


        
class graphene(lattice):
    def __init__(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                          {'tag': b'b', 'r0': [DX, DY]}]
        prim_vec = {'norm': 2 * DX, 'angle': 60}
        lattice.__init__(self, unit_cell=unit_cell, prim_vec=prim_vec)
        self.butterfly = np.array([])
        self.betas = np.array([])

    def triangle_zigzag(self, n):
        '''
        Triangular flake with zigzag terminations.

        :param: n. Int. Number of plackets along each edge. 
        '''
        check_n(n)
        self.get_lattice(n1=n, n2=n)
        ind = np.argwhere(self.coor['y'] < -2 * DX * self.coor['x'] + 3 * (n - 1)+1e-4)
        self.get_shape(ind)

    def hexagon_zigzag(self, n):
        '''
        Hexagonal flake with zigzag terminations.

        :param: n. Int. Number of plackets along each edge. 
        '''
        check_n(n)
        self.get_lattice(n1=2*n, n2=2*n)
        ind = np.argwhere((self.coor['y'] > -2 * DX  * self.coor['x'] + 3 * (n - 1) + 1e-4) &
                                   (self.coor['y'] < -2 * DX  * self.coor['x'] + 9 * n -1 - 1e-4))
        self.get_shape(ind)

    def triangle_armchair(self, n):
        '''
        Triangular flake with armchair terminations.

        :param: n. Int. Number of plackets along each edge. 
        '''
        check_n(n)
        self.get_lattice(n1=2*n, n2=2*n)
        ind = np.argwhere((self.coor['x'] > DX * (2*n-1)-1e-4) & 
                                   (self.coor['y'] >  0.5 / DX  * self.coor['x']  - n - 1e-4) &
                                   (self.coor['y'] < -0.5 / DX  * self.coor['x'] + 4 * n - 1 + 1e-4))
        self.get_shape(ind)

    def hexagon_armchair(self, n):
        '''
        Hexagonal flake with armchair terminations.

        :param: n. Int. Number of plackets along each edge. 
        '''
        check_n(n)
        nn = 3 * n - 2
        self.get_lattice(n1=2*nn, n2=2*nn)
        ind = np.argwhere((self.coor['x'] > DX * (2*nn-1)-1e-4) & 
                                   (self.coor['y'] >  0.5 / DX  * self.coor['x']  - nn - 1e-4) &
                                   (self.coor['y'] < -0.5 / DX  * self.coor['x'] + 4 * nn - 1 + 1e-4))
        self.get_shape(ind)
        self.coor['x'] -= self.coor['x'].min()
        ind = np.argwhere((self.coor['x'] < DX * 2*nn+0.1) & 
                                   (self.coor['y'] <  0.5 / DX  * self.coor['x']  +6*n - 4) &
                                   (self.coor['y'] >  -0.5 / DX  * self.coor['x']  +3*n-3))
        self.get_shape(ind)

    def square(self, n):
        '''
        Squared flake.

        :param: n. Int. Number of plackets along x. 
        '''
        check_n(n)
        n2 = int(1 + (2 * n * DX - 2) / 2.5)
        self.get_lattice(n1=n, n2=n2)
        self.sites = len(self.coor)
        ind = np.argwhere((self.coor['x'] > DX * (n-1)-0.1) & (self.coor['x'] < DX * 2*n-0.1))
        self.get_shape(ind)

    def circle(self, n):
        '''
        Circular flake.

        :param: n. Int. Number of plackets along x. 
        '''
        self.get_lattice(n1=n, n2=n)
        self.sites = len(self.coor)
        self.coor_center()
        ind = np.argwhere(self.coor['x'] ** 2 +  self.coor['y'] ** 2 < (0.75 * n) ** 2)
        self.get_shape(ind)
        
    def get_shape(self, ind):
        """
        Private method. Get the shape of the flake.
        """
        mask = np.zeros(self.sites, bool)
        mask[ind] = True
        self.coor = self.coor[mask]
        self.coor = np.sort(self.coor, order=('x', 'y'))   
        self.sites = len(self.coor)
        
class grapheneEig(system):
    def __init__(self, lat):
        system.__init__(self, lat)
        self.lat.coor['x'] -= np.round(self.lat.coor['x'].mean(), 2)
        self.lat.coor['y'] -= np.round(self.lat.coor['y'].mean(), 2)
        self.lat.coor = np.sort(self.lat.coor, order=('x', 'y'))

    def set_hop_linear_strain(self, t, beta):
        '''
        Set flake hoppings according to the linear trixial strain. 
        
        :param t: Hopping value without strain.
        :param beta: Strength of the strain.
        '''
        test_hop_linear_strain(t, beta)
        ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[1]-1e-4) &
                                   (self.vec_hop['d'] < self.dist_uni[1]+1e-4))
        ind_up = ind[ind[:, 1] > ind[:, 0]]  
        self.hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                                ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
        self.hop['n'] = 1
        self.hop['i'] = ind_up[:, 0]
        self.hop['j'] = ind_up[:, 1]
        self.hop['ang'] = self.vec_hop['a'][ind_up[:, 0], ind_up[:, 1]]
        self.hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], self.lat.coor['tag'][ind_up[:, 1]])
        # change angle (to get the correct strain)
        self.hop['ang'][self.hop['ang'] == 30] = -150
        x_center = .5 * (self.lat.coor['x'][ind_up[:, 0]] + self.lat.coor['x'][ind_up[:, 1]])
        y_center = .5 * (self.lat.coor['y'][ind_up[:, 0]] + self.lat.coor['y'][ind_up[:, 1]])
        self.hop['t'] = t * (1. + 0.25 * beta * (np.cos(PI / 180 * self.hop['ang']) * x_center +
                                                              np.sin(PI / 180 * self.hop['ang']) * y_center))
        # back to the former angle
        self.hop['ang'][self.hop['ang'] == -150] = 30

    def set_hop_nearest_hop(self, t, beta):
        '''
        Set flake hoppings according to the linear trixial strain. 
        
        :param t: Hopping value without strain.
        :param beta: Strength of the strain.
        '''
        test_hop_linear_strain(t, beta)
        ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[1]-1e-4) &
                                   (self.vec_hop['d'] < self.dist_uni[1]+1e-4))
        ind_up = ind[ind[:, 1] > ind[:, 0]]  
        self.hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                                ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
        self.hop['n'] = 1
        self.hop['i'] = ind_up[:, 0]
        self.hop['j'] = ind_up[:, 1]
        self.hop['ang'] = self.vec_hop['a'][ind_up[:, 0], ind_up[:, 1]]
        self.hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], self.lat.coor['tag'][ind_up[:, 1]])
        # change angle (to get the correct strain)
        
        # back to the former angle
        self.hop['ang'][self.hop['ang'] == -150] = 30

    def get_butterfly(self, t, N):
        ''''
        Get energies depending on strain.

        :param t: Hopping value.
        :param N: number of strain values between the max strains.
        '''
        test_hop_linear_strain(t, N)
        beta_lims = self.get_beta_lims()
        self.betas = np.linspace(beta_lims[0], beta_lims[1], N)
        self.butterfly = np.zeros((N, self.lat.sites))
        for i, beta in enumerate(self.betas):
            self.set_hop_linear_strain(t=1, beta=beta)
            self.get_ham(complex_transpose=True)
            self.butterfly[i] = LA.eigvalsh(self.ham.toarray())

    def get_beta_lims(self):
        '''
        Get the extremal values of strain keeping positive hoppings.
        '''
        beta_lims = np.zeros(2)
        yb_min_val = self.lat.coor['y'][self.lat.coor['tag'] == b'b'].min()
        yb_min = self.lat.coor['y'][self.lat.coor['y'] == yb_min_val][0]
        ym = 0.5*(2 * yb_min + 1)
        beta_lims[1] = -4./ym + 1e-6
        yb_max_val = self.lat.coor['y'][self.lat.coor['tag'] == b'a'].max()
        yb_max = self.lat.coor['y'][self.lat.coor['y'] == yb_max_val][0]
        ym = 0.5*(2 * yb_max - 1)
        beta_lims[0] = -4./ym + 1e-6
        print('Strain limits: {}'.format(beta_lims))
        return beta_lims

    def set_vortex(self):
        '''
        Set a vortex at the center of mass of the flake.
        '''
        ind = ((self.hop['ang'] == 90) & (self.lat.coor['y'][self.hop['i']] > -1.) & 
                (self.lat.coor['x'][self.hop['i']] <0) & (self.lat.coor['y'][self.hop['i']] < 0))
        self.hop['t'][ind] *= -1 
        if ind.sum() == 0:
            ind = ((self.hop['ang'] == 90) & (self.lat.coor['y'][self.hop['i']] > -2.) & 
                (self.lat.coor['x'][self.hop['i']] <0) & (self.lat.coor['y'][self.hop['i']] < 0))
            self.hop['t'][ind] *= -1 



lat = graphene()
lat.hexagon_armchair(n=5)
lat.plot()
plt.show()

