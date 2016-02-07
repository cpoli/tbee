from tbee.lattice import *
from tbee.plot import *
from tbee.system import *
from math import sqrt

PI = np.pi
ATOL = 1e-3
DX = 0.5 * sqrt(3)
DY = 0.5

 
#################################
# CLASS GRAPHENE
#################################

       
class grapheneLat(lattice):
    def __init__(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}, 
                          {'tag': b'b', 'r0': (DX, DY)}]
        prim_vec = [(2*DX, 0.), (DX, 1.5)]
        lattice.__init__(self, unit_cell=unit_cell, prim_vec=prim_vec)
        self.butterfly = np.array([])
        self.betas = np.array([])

    def triangle_zigzag(self, n):
        '''
        Triangular flake with zigzag terminations.

        :param: n. Int. Number of plackets along the edges. 
        '''
        error_handling.positive_int(n, 'n')
        self.get_lattice(n1=n+2, n2=n+2)
        self.boundary_line(cx=-sqrt(3), cy=-1, co=-3*n+2)

    def hexagon_zigzag(self, n):
        '''
        Hexagonal flake with zigzag terminations.

        :param: n. Int. Number of plackets along the edges. 
        '''
        error_handling.positive_int(n, 'n')
        self.get_lattice(n1=2*n, n2=2*n)
        self.boundary_line(cx=sqrt(3), cy=1, co=3*(n-1))
        self.boundary_line(cx=-sqrt(3), cy=-1, co=-9*n+2.5)

    def triangle_armchair(self, n):
        '''
        Triangular flake with armchair terminations.

        :param: n. Int. Number of plackets along the edges. 
        '''
        error_handling.positive_int(n, 'n')
        self.get_lattice(n1=2*n, n2=2*n)
        self.boundary_line(cx=1, cy=0, co=sqrt(3)/2*(2*n-1)-0.1)
        self.boundary_line(cx=-1/sqrt(3), cy=1, co=-n-0.1)
        self.boundary_line(cx=-1/sqrt(3), cy=-1, co=-4*n+0.1)

    def hexagon_armchair(self, n):
        '''
        Hexagonal flake with armchair terminations.

        :param: n. Int. Number of plackets along each edge. 
        '''
        error_handling.positive_int(n, 'n')
        nn = 3 * n - 2
        self.get_lattice(n1=2*nn, n2=2*nn)
        self.boundary_line(cx=1, cy=0, co=sqrt(3)/2* (2*nn-1)-.1)
        self.boundary_line(cx=-1/sqrt(3), cy=1, co=-nn -.1)
        self.boundary_line(cx=-1/sqrt(3), cy=-1, co=-4*nn +1-.1)
        self.coor['x'] -= self.coor['x'].min()
        self.boundary_line(cx=-1, cy=0, co=-DX * 2*nn-0.1)
        self.boundary_line(cx=1/sqrt(3), cy=1, co=3*n -3.)
        self.boundary_line(cx=1/sqrt(3), cy=-1, co=-6*n+4)

    def square(self, n):
        '''
        Squared flake.

        :param: n. Int. Number of plackets along x. 
        '''
        error_handling.positive_int(n, 'n')
        n2 = int(1.5*DX*n)
        self.get_lattice(n1=2*n, n2=n2)
        self.boundary_line(cx=1, cy=0, co=DX*(2*n-2)+0.1)
        self.boundary_line(cx=-1, cy=0, co=-DX*(4*n)+0.5)

    def circle(self, n):
        '''
        Circular flake.

        :param: n. Int. Number of plackets along the diameter. 
        '''
        error_handling.positive_int(n, 'n')
        self.get_lattice(n1=2*n, n2=2*n)
        self.sites = len(self.coor)
        self.center()
        if n % 2 == 0:
            self.shift_x(shift=-DX)
        self.ellipse_in(rx=DX*(n+1), ry=DX*(n+1), x0=0., y0=0.)
        self.remove_dangling()
        
class grapheneSys(system):
    def __init__(self, lat):
        system.__init__(self, lat)

    def set_hop_linear_strain(self, t, beta):
        '''
        Set nearest neighbors hoppings according to the linear trixial strain. 
        
        :param t: Hopping value without strain.
        :param beta: Strength of the strain.
        '''
        error_handling.number(t, 't')
        error_handling.real_number(beta, 'beta')
        self.get_distances()
        ind = np.argwhere(np.isclose(self.dist_uni[1], self.vec_hop['dis'], atol=ATOL))
        ind_up = ind[ind[:, 1] > ind[:, 0]]  
        self.hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                                        ('ang', 'f8'), ('tag', 'S2'), ('t', 'c16')])
        self.hop['n'] = 1
        self.hop['i'] = ind_up[:, 0]
        self.hop['j'] = ind_up[:, 1]
        self.hop['ang'] = self.vec_hop['ang'][ind_up[:, 0], ind_up[:, 1]]
        # change angle (to get the correct strain)
        self.hop['ang'][np.isclose(30., self.hop['ang'], ATOL)] = -150.
        self.hop['ang'][np.isclose(150., self.hop['ang'], ATOL)] = - 30.
        x_center = .5 * (self.lat.coor['x'][ind_up[:, 0]] + self.lat.coor['x'][ind_up[:, 1]])
        y_center = .5 * (self.lat.coor['y'][ind_up[:, 0]] + self.lat.coor['y'][ind_up[:, 1]])
        self.hop['t'] = t * (1. + 0.25 * beta * (np.cos(PI / 180 * self.hop['ang']) * x_center +
                                                                np.sin(PI / 180 * self.hop['ang']) * y_center))
        # back to the former angle
        self.hop['ang'][np.isclose(-150., self.hop['ang'])] = 30.
        self.hop['ang'][np.isclose(-30., self.hop['ang'])] = 150.

    def get_butterfly(self, t, N):
        ''''
        Get energies depending on strain.

        :param t: Unstrained hopping value.
        :param N: number of strain values between min and max strains.
        '''
        error_handling.number(t, 't')
        error_handling.positive_int(N, 'N')
        beta_lims = self.get_beta_lims()
        self.betas = np.linspace(beta_lims[0], beta_lims[1], N)
        self.butterfly = np.zeros((N, self.lat.sites))
        for i, beta in enumerate(self.betas):
            self.set_hop_linear_strain(t=1, beta=beta)
            self.get_ham()
            self.butterfly[i] = LA.eigvalsh(self.ham.toarray())

    def get_beta_lims(self):
        '''
        Get the extremal values of strain keeping positive hoppings.
        '''
        beta_lims = np.zeros(2)
        yb_min_val = self.lat.coor['y'][self.lat.coor['tag'] == b'b'].min()
        yb_min = self.lat.coor['y'][self.lat.coor['y'] == yb_min_val][0]
        ym = 0.5 * (2 * yb_min + 1)
        beta_lims[1] = -4. / ym + 1e-6
        yb_max_val = self.lat.coor['y'][self.lat.coor['tag'] == b'a'].max()
        yb_max = self.lat.coor['y'][self.lat.coor['y'] == yb_max_val][0]
        ym = 0.5 * (2 * yb_max - 1)
        beta_lims[0] = -4. / ym + 1e-6
        print('Strain limits: {}'.format(beta_lims))
        return beta_lims
