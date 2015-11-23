"""
DOUBLE CHECK DANGLING
"""
import numpy as np
import matplotlib.pyplot as plt


PI = np.pi


def test_ri(ri):
    '''
    Check parameter *ri*.

    :raises TypeError: Parameters *ri* must be a list.
    :raises ValueError: Parameter *ri* must be a container of 2D coordinates.
      Example: ri = [[0., 1.], [2., 3.]]  -> 2 points in 2D.
    :raises ValueError: Parameter *ri* must contain real numbers.
    '''
    error_not_2d = '\n\nParameter ri must be a container of 2D coordinates.\
                             Example: ri = [[0., 1.], [2., 3.]] -> 2 points.\n.'
    if not isinstance(ri, list):
        raise TypeError('\n\nParameters ri must be a list.\n')
    if ri == [] or not all([isinstance(r, list) and len(r) == 2 for r in ri]):
        raise ValueError(error_not_2d)
    if not all([isinstance(x, (int, float)) and isinstance(y, (int, float)) for [x, y] in ri]):
        raise ValueError(error_not_2d)


def test_tags(tags):
    '''
    Check parameter *tags*.

    :raises TypeError: Parameter *tags* must be a list.
    :raises ValueError: Parameter *tags* must be a list of binary chars.
    '''
    if not isinstance(tags, list):
        raise TypeError('\n\nParameter tags must be a list.\n')
    if not all([isinstance(tag, bytes) and len(tag) == 1 for tag in tags]):
        raise ValueError('\n\nParameter tags must be a list of binary chars.\n')


def test_nor(nor):
    '''
    Check parameter *nor*.

    :raises TypeError: Parameter nor must be a positive number.
    '''
    if not isinstance(nor, (int, float)) or nor < 0:
        raise TypeError('\n\nParameter nor is not a positive number.\n')


def test_ang(ang):
    '''
    Check parameter *ang*.

    :raises TypeError: Parameter ang  must be a number.
    '''
    if not isinstance(ang, (int, float)):
        raise TypeError('\n\nParameter ang must be a real number.\n')

def test_remove_sites(index, sites):
    '''
    Check method *remove_sites*.

    :raises TypeError: Parameter index must be a list.
    :raises IndexError: Parameter index must be a list
      of integers between 0 and sites
    '''
    if not isinstance(index, list):
        raise TypeError('\n\nParameter index must be a list.\n')
    print(index)
    if not all(isinstance(i, int) for i in index):
        raise IndexError('\n\nParameter index must be a list\
                                 of integers between 0 and sites - 1.\n')
    if not all(-1 < i < sites for i in index):
        raise IndexError('\n\nParameter index must be a list\
                                    of integers between 0 and sites - 1.\n')


def test_plt_lattice(tags, ms, fs, colors, plt_index, figsize):
    '''
    Check method *remove_sites*.

    :raises TypeError: Parameter *ms* must be a positive integer.
    :raises TypeError: Parameter *fs* must be a positive integer.
    :raises TypeError: Parameter *plt_index* must be a bool.
    :raises TypeError: Parameter *figsize* must be list/tuple.
    :raises ValueError: Parameter *figsize* must contain positive numbers.
    :raises ValueError: Parameter *colors* must contain 
      RGB or matplotlib colors.
    '''
    if not isinstance(ms, int) or ms < 1:
        raise TypeError('\n\nParameter ms must be a positive integer.\n')
    if not isinstance(fs, int) or fs < 1:
        raise TypeError('\n\nParameter fs must be a positive integer.\n')
    if not isinstance(plt_index, bool):
        raise TypeError('\n\nParameter plt_index must be a bool.\n')
    if figsize is not None:
        if not isinstance(figsize, (list, tuple)):
            raise TypeError('\n\nParameter figsize must be a list/tuple\n')
        if len(figsize) != 2:
            raise ValueError('\n\nParameter figsize must be of length 2\n')
        if isinstance(figsize, (list, tuple)):
            if any([val < 0 for val in figsize]):
                raise ValueError('\n\nParameter figsize must contain\
                                        positive numbers.\n')
    if not isinstance(colors, list):
        colors = ['b', 'r', 'g', 'y', 'm', 'k', 'c']
    else:
        if len(colors) < len(tags):
            raise ValueError('\n\nlength of colors must be a least equal to\
                                        length of tags.\n')
        plt_colors = ['b', 'r', 'g', 'y', 'm', 'k', 'c']
        for c in colors:
            if not (c[0] == '#' and len(c) == 7) and c not in plt_colors:
                raise ValueError('\n\nParameter colors must be a list/tuple of\
                                          of hexagonal colors or matplotlib colors:\
                                          ["b", "g", "r", "c", "m", "y", "k"]\n')
    return colors, figsize

class latticeTB(object):
    '''
    Build up 1D or 2D lattice.
    '''
    def __init__(self, ri, tags, nor, ang):
        '''
        :param tags: The tags (binary chars) associated to each sublattice.
        :param ri: 2D Coordinates. Initial positions of sites
          belonging to the unit cell.
        :param nor: Primitive vectors norm.
        :param ang: Angle between both primitive vectors.
        '''
        test_ri(ri)
        test_tags(tags)
        if len(ri) != len(tags):
            raise ValueError('\n\nParameters ri and tags must have the same length\n')
        test_ang(ang)
        test_nor(nor)
        self.tags = np.array(tags)  # name of the tags
        self.ri = np.array(ri)  # positions of the sites within the unit cell
        self.nx, self.ny = 0, 0  # no of unit cells along x and y
        self.ri = ri  # initial positions of the sites
        self.nor = nor  # primitive vectors norm
        self.ang = ang  # angle between primitive vectors
        self.coor = np.array([], dtype=[('x', 'f16'), ('y',  'f16') , ('tag', 'S1')])
        self.sites = 0  # Site number

    def test_coor(self):
        '''
        Check if *get_lattice* has been called (*self.coor* not empty).
        :raises RuntimeError: *self.coor* empty.
        '''
        if self.coor.size == 0:
            raise RuntimeError('\n\nRun method get_lattice first\n')

    def get_lattice(self, nx, ny):
        '''
        Get the lattice positions

        :param nx: Number of unit cells along :math:`x`.
        :param ny: Number of unit cells along :math:`y`.
        '''
        if not isinstance(nx, int) or nx < 0:
            raise ValueError('\n\nParameter nx is not a positive integer.\n')
        if not isinstance(ny, int) or ny < 0:
            raise ValueError('\n\nParameter ny is not a positive integer.\n')

        self.coor = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        self.nx, self.ny = nx, ny
        x = self.nor * np.arange(nx)
        y = np.sin(self.ang) * self.nor * np.arange(ny)
        xx, yy = np.meshgrid(x, y)
        if self.ang and self.ang != PI/2:
            xx += yy
        elif self.ang:
            xx += yy / np.tan(self.ang)
        xx = np.ravel(xx)
        yy = np.ravel(yy)
        coor_tag = np.zeros(len(xx), dtype=[('x', 'f16'), ('y',  'f16') , ('tag', 'S1')])
        for i, t in enumerate(self.tags):
            coor_tag['x'] = xx + self.ri[i][0]
            coor_tag['y'] = yy + self.ri[i][1]
            coor_tag['tag'] = t
            self.coor = np.concatenate([self.coor, coor_tag])
        self.coor = np.sort(self.coor, order=('y', 'x'))        
        self.coor['x'] = self.coor['x'].round(6)
        self.coor['y'] = self.coor['y'].round(6)
        self.sites = len(self.coor['tag'])

    def remove_sites(self, index):
        '''
        Remove sites define their indices

        :param index: Array of site indices to be removed.
        '''
        self.test_coor()
        test_remove_sites(index, self.sites)
        mask = np.ones(self.sites, bool)
        mask[index] = False
        self.coor = self.coor[mask]
        self.sites = self.coor.size

    def remove_dangling(self, nor_bond):
        '''
        Remove dangling sites
        (sites connected with just another site).

        :param nor_bond: Norm of the bonds to be checked.
        '''
        self.test_coor()
        if not isinstance(nor_bond, (int, float)) or nor_bond < 0:
            raise TypeError('\n\nParameter nor_bond is not a positive number.\n')

        while True:
            dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
            dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
            dis = np.sqrt(dif_x ** 2 + dif_y ** 2)
            ind = np.argwhere((dis > nor_bond-1e-9) & (dis < nor_bond+1e-9))
            dang = []
            for i in range(self.sites):
                if (ind[:, 0] == i).sum() == 1:
                    dang.append(i)
            self.coor = np.delete(self.coor, dang, axis=0)
            self.sites -= len(dang)
            if dang == []:
                break

    def plt_lattice(self, ms=30, fs=20, colors=None, plt_index=False, figsize=None):
        '''
        Plot lattice.

        :param colors: Sublattice colors.
        :param ms: Default value 30. Markersize.
        :param lw: Default value 5. Linewidth of the hoppings.
        :param plt_index: Default value False. Plot site indices.
        :param figsize: Default value None. Figsize.

        :returns:
            * **fig** -- Figure.
        '''
        self.test_coor()
        colors, figsize = test_plt_lattice(self.tags, ms, fs, colors, plt_index, figsize)
        fig, ax = plt.subplots(figsize=figsize)
        # plot sites
        for c, t in zip(colors, self.tags):
            str_tag = t.upper()   # .decode('ascii')
            plt.plot(self.coor['x'][self.coor['tag'] == t],
                        self.coor['y'][self.coor['tag'] == t], 'o', color=c,
                        ms=ms, label=str_tag, markeredgecolor='none')
        ax.axis('off')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(self.coor['x'])-0.5, np.max(self.coor['x'])+0.5])
        ax.set_ylim([np.min(self.coor['y'])-0.5, np.max(self.coor['y'])+0.5])
        # plot indices
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sites)]
            for l, x, y in zip(indices, self.coor['x'], self.coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        return fig
