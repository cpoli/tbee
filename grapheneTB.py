'''
    TO BE MODIFIED

    **grapheneTB** module models graphene flakes.
    Flakes can be submitted to triaxial linear strain.
    
    **grapheneTB** is available at https://github.com/cpoli/TB
'''

from latticeTB import *
from plotTB import *
from eigTB import *
from propagationTB import *
from math import sqrt, pi

dx, dy = sqrt(3)/2., 0.5


def rec(nx, ny):
    '''
    Graphene rectangular flake.
    Used to build up the other flakes. 

      :param nx: Number of plaquettes along :math:`x`. 
      :param ny: Number of plaquettes along :math:`x`. 
    '''
    sites_x = 2 * nx + 3
    coor = np.zeros(2*sites_x*ny, dtype={'names': ['x', 'y', 'tag'], 
                                                                'formats': ['f8', 'f8', 'S1']})
    row = dx*np.arange(0, sites_x)
    sites, y = 0, 0.
    for i in range(2*ny):
        coor['x'][sites: sites+sites_x] = row
        if i % 2 == 0:
            coor['y'][0+sites: sites+sites_x: 2] = y
            coor['y'][1+sites: sites+sites_x: 2] = y + dy
            coor['tag'][0+sites: sites+sites_x: 2] = 'a'
            coor['tag'][1+sites: sites+sites_x: 2] = 'b'
        else:
            coor['y'][0+sites: sites+sites_x: 2] = y + dy
            coor['y'][1+sites: sites+sites_x: 2] = y
            coor['tag'][0+sites: sites+sites_x: 2] = 'b'
            coor['tag'][1+sites: sites+sites_x: 2] = 'a'
        y += 1.5
        sites += sites_x
    return coor, sites


def square(n):
    '''
    Graphene squared flake.

    :param n: Number of plaquettes along :math:`x`. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx, ny = n, int(((n+1)*sqrt(3)+1)/3)
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    ind = []
    coor = remove(ind, coor_rec, sites_rec)
    coor_rec['x'] -= np.mean(coor_rec['x'])
    coor_rec['y'] -= np.mean(coor_rec['y'])
    return coor_rec, coor_rec, ind, nx, ny


def circle(n):
    '''
    Graphene circular flake.

    :param n: Number of plaquettes along radius. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx = 2*n
    ny = int(1.3*n)
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    coor_rec['x'] -= np.mean(coor_rec['x'])
    coor_rec['y'] -= np.mean(coor_rec['y'])
    ind = np.argwhere(np.sqrt(coor_rec['x']**2+coor_rec['y']**2) > 2*dx*n)
    coor = remove(ind, coor_rec, sites_rec)
    return coor, coor_rec, ind, nx, ny


def tri_zz(n):
    '''
    Graphene triangular flake with zigzag terminations.

    :param n: Number of plaquettes along one edge. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx = n
    if nx % 2 == 0:
        ny = nx//2 + 1
    else:
        ny = nx//2 + 2
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    coor_rec['x'] -= np.mean(coor_rec['x'])
    L = sqrt(3)*(n+1)
    h = 1.5*(n+1)
    l = sqrt(h**2+0.25*L**2)
    ind = np.argwhere(coor_rec['y'] > h-2*np.abs(coor_rec['x'])*h/L+0.01)
    coor_rec['y'] -= np.mean(coor_rec['y'])
    coor = remove(ind, coor_rec, sites_rec)
    return coor, coor_rec, ind, nx, ny


def hex_zz(n):
    '''
    Graphene hexagonal flake with zigzag terminations.

    :param n: Number of plaquettes along one edge. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx = 2*n
    ny = n + 1
    if nx %2 == 0:
        nx -= 1
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    Lx = sqrt(3) * n
    Ly = 1.5*(n-1) + 1
    coor_rec['x'] -= np.mean(coor_rec['x'])
    coor_rec['y'] -= np.mean(coor_rec['y'])
    ind = np.argwhere((np.abs(coor_rec['y']) > Ly) |
                      (np.abs(coor_rec['y']) > -2*Ly/Lx*np.abs(coor_rec['x'])+2*Ly))
    coor = remove(ind, coor_rec, sites_rec)
    return coor, coor_rec, ind, nx, ny


def tri_ac(n):
    '''
    Graphene triangular flake with armchair terminations.

    :param n: Number of plaquettes along :math:`y`. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx = n
    if nx %2 == 0:
        ny=nx//2+2
    else:
        ny=nx//2+3
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    L = sqrt(3)*(n+1)
    h = 1.5*(n+1)
    l = sqrt(h**2+0.25*L**2)
    coor_rec['y'] -= np.mean(coor_rec['y'])
    ind = np.argwhere(coor_rec['x'] > h-2*np.abs(coor_rec['y'])*h/L+0.01)
    coor_rec['y'] -= np.mean(coor_rec['y'])
    coor_rec['x'] -= np.mean(coor_rec['x'])
    coor = remove(ind, coor_rec, sites_rec)
    coor['y'] -= np.mean(coor['y'])
    coor['x'] -= np.mean(coor['x'])
    return coor, coor_rec, ind, nx, ny


def hex_ac(n):
    '''
    Graphene hexagonal flake with armchair terminations.

    :param n: Number of plaquettes along :math:`y`. 

    :returns:
      * **coor** -- Flake coordinates. 
      * **coor_rec** -- Rectangular flake coordinates. 
      * **ind** -- Indices to be removed. 
      * **nx** -- Rectangular flake, number of plaquettes along :math:`x`. 
      * **ny** -- Rectangular flake, number of plaquettes along :math:`y`. 
    '''
    nx = 2*n
    ny = int(1.2*n)+1
    if nx % 2 == 0:
        nx -= 1
    coor_rec, sites_rec = rec(nx=nx, ny=ny)
    Lx = sqrt(3)*n
    Ly = 1.5*(n-1)+1
    coor_rec['x'] -= np.mean(coor_rec['x'])
    coor_rec['y'] -= np.mean(coor_rec['y'])
    ind = np.argwhere((np.abs(coor_rec['x']) > Ly) |
                      (np.abs(coor_rec['x']) > -2*Ly/Lx*np.abs(coor_rec['y'])+2*Ly))
    coor = remove(ind, coor_rec, sites_rec)
    return coor, coor_rec, ind, nx, ny


def remove(ind, coor_rec, sites_rec):
    '''
    Remove sites from the rectangular flake to get the desired flake.

    :param ind: Indices to be removed. 
    :param coor_rec: Rectangular flake coordinates. 
    :param sites_rec: Rectangular flake number of sites. 

    :returns:
      * **coor** -- Flake coordinates. 
    '''
    mask = np.ones(sites_rec, bool)
    mask[ind] = False
    coor = coor_rec[mask]
    coor['x'] -= np.mean(coor['x'])
    coor['y'] -= np.mean(coor['y'])
    return coor


class grapheneTB():
    '''
    Build up a graphene flake.
    '''
    def __init__(self):
        self.tags = np.array([b'a', b'b'])
        self.sites = -1  # Sites flake.
        self.sites_rec = -1  # Sites rectangular flake.
        self.coor = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.coor_rec = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.ind = np.array([])  # indices to build up the flake
        self.nx = -1  # Rectangular flake. number of plaquettes along :math:`x`
        self.ny = -1  # Rectangular flake. number of plaquettes along :math:`y`
        self.ind_del = -1  # Site indices to be deleted. 

    def get_lattice(self, flake, n):
        '''
        Get graphene flake.

        :param flake: Flake geometry.
        '''
        self.coor, self.coor_rec, self.ind, self.nx, self.ny = flake(n)
        self.sites = len(self.coor['tag'])
        self.sites_rec = len(self.coor_rec['tag'])
        print('sites', self.sites)

class plotGraphene(plotTB):
    '''
    Child of the class **plotTB**. Dedicated to Graphene. 
    '''
    def __init__(self, coor, tags, coor_hop=[]):
        plotTB.__init__(self, coor, tags)
        self.coor_hop = np.array([], {'names': ['x', 'y', 'tag'], 
                                                     'formats': ['f8', 'f8', 'S1']})

    def plt_butterfly(self, butterfly, betas, en_lims=[-2., 2.], fs=20):
        '''
        Plot energies depending on strain.
        '''
        if not butterfly.any():
            raise Exception('\n\nRun method get_butterfly() first.\n')
        i_beta_min = np.argmin(np.abs(betas))
        ind_en = np.argwhere((butterfly[:, i_beta_min] > en_lims[0]) & 
                                              (butterfly[:, i_beta_min] < en_lims[1]))
        ind_en = np.ravel(ind_en)
        fig, ax = plt.subplots()
        plt.title('Energies depending on strain', fontsize=fs)
        plt.xlabel(r'$\beta/\beta_{max}$', fontsize=fs)
        plt.ylabel('$E$', fontsize=fs)
        plt.xticks((-betas[-1], -.5*betas[-1], 0, .5*betas[-1], betas[-1]), 
                        fontsize=fs)
        plt.yticks(np.arange(en_lims[0], en_lims[1]+1e-2), fontsize=fs)
        ax.set_xticklabels(('-1', '-1/2', '0', '1/2', '1'))
        plt.xlim([betas[0], betas[-1]])
        plt.ylim(en_lims)
        print(butterfly.shape)
        print(butterfly[ind_en, 0])
        no_en = len(ind_en)
        for i in ind_en:
            plt.plot(betas, butterfly[i, :], 'b')
        fig.set_tight_layout(True)
        plt.draw()
        return fig
        
    def plt_lattice(self, coor, hop=[], ms=30, lw=5, fs=20, c=3., label=None, figsize=None):
        '''
        :param coor: Coordinates.
        :param hop: Default value []. Hoppings given by the class **eigTB**. 
            If not empty, plot the hoppings.
        :param ms: Default value 30. Markersize. 
        :param lw: Default value 5. Linewidth of the hoppings.
        :param fs: Default value 20. Fontsize.
        :param c: Default value 3. Coefficient, Bonds linewidth given by c*hop.
        :param label: Default value False. Plot site labels.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** Figure.
        '''
        fig, ax = plt.subplots(figsize=figsize)
        # bonds
        if hop != []:
            for i in range(len(hop)): 
                plt.plot([coor['x'][hop['i'][i]], coor['x'][hop['j'][i]]],
                                [coor['y'][hop['i'][i]], coor['y'][hop['j'][i]]],
                                'k', lw=c*hop['t'][i].real)
        # sites
        for c, t in zip(self.colors, self.tags):
            plt.plot(coor['x'][coor['tag'] == t],
                    coor['y'][coor['tag'] == t],
                    'o', color=c, ms=ms, markeredgecolor='none')
        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_xlim([np.min(coor['x'])-0.5, 
                      np.max(coor['x'])+0.5])
        ax.set_ylim([np.min(coor['y'])-0.5, 
                      np.max(coor['y'])+0.5])
        # labels
        if label:
            labels = ['{}'.format(i) for i in range(len(coor['tag']))]
            for l, x, y in zip(labels, coor['x'], coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        plt.draw()
        return fig

    def plt_dispersion(self, t, N=200, fs=20):
        '''
        Plot dispersion relation of Graphene.

        :param t: Hopping value.
        :param N: Default value 100. Number of points along each direction.
        :param fs: Default value 20. Fontsize.
        '''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        Ek = np.zeros((N, N, 2), 'f16')
        kx = np.linspace(-pi, pi, N)
        ky = np.linspace(-pi, pi, N)
        kx, ky = np.meshgrid(kx, ky) 
        en = t * np.sqrt(3. + 2. * np.cos(2* dx * ky) + 4. * np.cos(dx * ky) * 
                                np.cos(np.sqrt(3) * dx * kx))
        Ek[:, :, 0] = - np.sqrt(en)
        Ek[:, :, 1] = np.sqrt(en)
        alpha = 0.8
        ax.plot_surface(kx, ky, Ek[:, :, 0], linewidth=0., antialiased=True, 
                shade=True, alpha=alpha, color='b')
        ax.plot_surface(kx, ky, Ek[:, :, 1], linewidth=0., antialiased=True, 
                shade=True, alpha=alpha, color='r')
        plt.title('$E(\mathbf{k})$', fontsize=fs)
        ax.set_xlabel('$k_x$', fontsize=fs)
        ax.set_ylabel('$k_y$', fontsize=fs)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([0])
        plt.draw()


class eigGraphene(eigTB):
    '''
    Child of the class **eigTB**. Dedicated to Graphene flakes. 
    '''
    def __init__(self, coor, coor_rec, tags):
        eigTB.__init__(self, coor, tags)
        self.beta = 0  # strain strength
        self.beta_lims = self.get_beta_lims()  # strain limits keeping the positive hops 
        self.coor_rec = coor_rec  # Rectangular flake coordinates
        self.hop_rec = np.array([], {'names': ['i', 'j', 't'], 'formats': ['u4', 'u4', 'c16']})
        self.coor_hop = np.array([], {'names': ['i', 'j', 't'], 'formats': ['u4', 'u4', 'c16']})
        self. beta = 0.  # strain strength
        self.betas = np.array([])  # array of betas to get the butterfly
        self.butterfly = np.array([])  # store the energies of the butterfly
        self.ind_del = np.array([])  # store the site indices to be deleted

    def get_hop(self, t1, t2, t3, coor):
        '''
        Do not call this method: call *get_hops*.
        '''
        sites = len(coor['tag'])
        pos_y = coor['y']
        dif_y = pos_y.reshape(sites, 1) - pos_y
        pos_x = coor['x']
        dif_x = pos_x.reshape(sites, 1) - pos_x
        ind1 = np.argwhere((dif_y == 1.) & (dif_x == 0.))
        ind2 = np.argwhere((dif_y == -0.5) & (dif_x > -dx-0.1) & (dif_x < -dx+0.1))
        ind3 = np.argwhere((dif_y == -0.5) & (dif_x > dx-0.1) & (dif_x < dx+0.1))
        bonds_y, bonds_x = len(ind1), len(ind2)
        bonds = bonds_y + 2*bonds_x
        hop = np.zeros(bonds, dtype={'names': ['i', 'j', 't', 'tag'],
                                                         'formats': ['i4', 'i4', 'f8', 'S1']})
        # t1
        hop['i'][0: bonds_y] = ind1[:, 0]
        hop['j'][0: bonds_y] = ind1[:, 1]
        hop['tag'][0: bonds_y] = 'u'
        hop['t'][0: bonds_y] = t1 * \
          (1. + 0.5 * dy * self.beta * (coor['y'][ind1[:, 0]] + coor['y'][ind1[:, 1]]))
        # t2
        hop['i'][bonds_y: bonds_y + bonds_x] = ind2[:, 0]
        hop['j'][bonds_y: bonds_y + bonds_x] = ind2[:, 1]
        hop['tag'][bonds_y: bonds_y + bonds_x] = 'v'
        hop['t'][bonds_y: bonds_y + bonds_x] = t2 * \
          (1. + 0.25 * self.beta * (- dx * (coor['x'][ind2[:, 0]] + coor['x'][ind2[:, 1]])
                                               - dy * (coor['y'][ind2[:, 0]] + coor['y'][ind2[:, 1]])))
        # t3
        hop['i'][bonds_y + bonds_x: bonds] = ind3[:, 0]
        hop['j'][bonds_y + bonds_x: bonds] = ind3[:, 1]
        hop['tag'][bonds_y + bonds_x: bonds] = 'w'
        hop['t'][bonds_y + bonds_x: bonds] = t3 * \
          (1. + 0.25 * self.beta * (+ dx * (coor['x'][ind3[:, 0]] + coor['x'][ind3[:, 1]])
                                                - dy * (coor['y'][ind3[:, 0]] + coor['y'][ind3[:, 1]])))
        return hop

    def get_hops(self, t1, t2, t3, beta=0., get_hop_rec=False):
        '''
        Get hoppings of flake and rectangular flake.
        Needed only if the lattice will be plotted in hopping space.

        :param t1: Hopping value.
        :param t2: Hopping value. 
        :param t3: Hopping value. 
        :param beta: Default value 0. Strain strength. 
        :param get_hop_rec: Default value False.
        '''
        self.beta = beta
        if get_hop_rec:
            self.hop_rec = self.get_hop(t1, t2, t3, self.coor_rec)
        self.hop = self.get_hop(t1, t2, t3, self.coor)

    def get_coor_hop(self, nx, ny, coor_rec, ind):
        '''
        Get coordinates in hopping space.

        :param coor_rec: Rectangular flake coordinates.
        :param nx: Rectangular flake, number of plaquettes along :math:`x`. 
        :param ny: Rectangular flake, number of plaquettes along :math:`y`. 
         :param ind: Indices to be removed. 
        '''
        sites, bonds = (2*nx+3)*2*ny, (2*nx+2)*2*ny
        sites_x, bonds_x = 2*nx+3, 2*nx+2
        coor_hop = np.zeros(sites, dtype={'names': ['x', 'y', 'tag'],
                                                                'formats': ['f4', 'f4', 'S1']})
        pos_y = coor_rec['y']
        dif_y = pos_y.reshape(sites, 1) - pos_y
        pos_x = coor_rec['x']
        dif_x = pos_x.reshape(sites, 1) - pos_x
        bonds1 = ((dif_y == 1.) & (dif_x == 0.)).sum()
        bonds2 = ((dif_y == -0.5) & (dif_x > -dx-0.1) & (dif_x < -dx+0.1)).sum()
        bonds3 = ((dif_y == -0.5) & (dif_x > dx-0.1) & (dif_x < dx+0.1)).sum()
        hop = self.hop_rec
        bonds_v, bonds_w = bonds1, bonds1+bonds2
        ty = np.unique(self.hop_rec[0: bonds1]['t'])
        if len(ty) == 1:
            ty = ty*np.ones(2*ny-1)
        t = np.zeros(bonds_x, dtype={'names': ['x', 'y'], 'formats': ['f8', 'f8']})
        sites = 0
        k = 0
        for i in range(2*ny):
            if i % 2 == 0:
                t['x'][0::2] = dx * hop['t'][bonds_v: bonds_v+bonds_x//2]
                t['y'][0::2] = dy * hop['t'][bonds_v: bonds_v+bonds_x//2]
                t['x'][1::2] = dx * hop['t'][bonds_w: bonds_w+bonds_x//2]
                t['y'][1::2] = - dy * hop['t'][bonds_w: bonds_w+bonds_x//2]
                coor_hop['tag'][0+sites: sites+sites_x:2] = 'a'
                coor_hop['tag'][1+sites: sites+sites_x:2] = 'b'
            else:
                t['x'][1::2] = dx * hop['t'][bonds_v: bonds_v+bonds_x//2]
                t['y'][1::2] = dy * hop['t'][bonds_v: bonds_v+bonds_x//2]
                t['x'][0::2] = dx * hop['t'][bonds_w: bonds_w+bonds_x//2]
                t['y'][0::2] = - dy * hop['t'][bonds_w: bonds_w+bonds_x//2]
                coor_hop['tag'][1+sites: sites+sites_x: 2] = 'a'
                coor_hop['tag'][0+sites: sites+sites_x: 2] = 'b'
            coor_hop['x'][sites+1: sites+sites_x] = coor_hop['x'][sites] \
                + np.cumsum(t['x'])
            coor_hop['y'][sites+1: sites+sites_x] = coor_hop['y'][sites] \
                + np.cumsum(t['y'])
            if i < 2*ny-1:
                if i % 2 == 0:
                    coor_hop['x'][sites+sites_x] = coor_hop['x'][sites+1] \
                        - dx * hop['t'][bonds_w+bonds_x//2]
                    coor_hop['y'][sites+sites_x] = coor_hop['y'][sites+1] \
                        + ty[k] + dy * hop['t'][bonds_w+bonds_x//2]
                    coor_hop['tag'][sites+sites_x] = 'b'
                else:
                    coor_hop['x'][sites+sites_x] = coor_hop['x'][sites]
                    coor_hop['y'][sites+sites_x] = coor_hop['y'][sites]+ty[k]
                    coor_hop['tag'][sites+sites_x] = 'a'
            k += 1
            sites += sites_x
            bonds_v += nx + 1
            bonds_w += nx + 1
        coor_hop = remove(ind, coor_hop, sites)
        coor_hop['x'] -= coor_hop['x'].mean()
        coor_hop['y'] -= coor_hop['y'].mean()
        coor_hop = np.delete(coor_hop, self.ind_del) 
        self.coor_hop = coor_hop

    def get_beta_lims(self):
        '''
        Get the extremal values of strain keeping positive hoppings.
        '''
        pos_y = self.coor['y']
        dif_y = pos_y.reshape(self.sites, 1) - pos_y
        pos_x = self.coor['x']
        dif_x = pos_x.reshape(self.sites, 1) - pos_x
        ind1 = np.argwhere((dif_y == 1.) & (dif_x == 0.))
        beta_lims = np.zeros(2)
        ym = 0.5*(self.coor['y'][ind1[-1][0]] + self.coor['y'][ind1[-1][1]])
        beta_lims[0] = -2./ym + 1e-2
        ym = 0.5*(self.coor['y'][ind1[0][0]] + self.coor['y'][ind1[0][1]])
        beta_lims[1] = -2./ym - 1e-2
        return beta_lims

    def get_butterfly(self, t1, t2, t3, N):
        ''''
        Get energies depending on strain.

        :param t1: Hopping value.
        :param t2: Hopping value. 
        :param t3: Hopping value. 
        :param N: Strain value number.
        '''
        import scipy.sparse.linalg as sla

        self.get_beta_lims()
        self.betas = np.linspace(self.beta_lims[0], self.beta_lims[1], N)
        self.butterfly = np.zeros((self.sites, N))
        for i, beta in enumerate(self.betas):
            self.get_hops(t1, t2, t3, beta, get_hop_rec=False)
            self.get_ham()
            self.butterfly[:, i] = LA.eigvalsh(self.ham.toarray())

    def del_sites(self, ind_del):
        ''''
        Delete sites.

        :param ind_del: Site indices to be deleted
        '''
        self.ind_del = np.array([ind_del])
        self.coor = np.delete(self.coor, self.ind_del) 
        if self.coor_hop != []:
            self.coor_hop = np.delete(self.coor_hop, self.ind_del) 
        h = self.ham.toarray()
        h = np.delete(h, self.ind_del, axis=0)
        h = np.delete(h, self.ind_del, axis=1)
        self.ham = sparse.csr_matrix(h)
        self.sites -= len(self.ind_del)

