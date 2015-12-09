import numpy as np
import matplotlib.pyplot as plt


PI = np.pi


def test_unit_cell(unit_cell):
    '''
    Check parameter *unit_cell*.

    :raises TypeError: Parameter unit_cell must be a list.
    :raises KeyError: Dictionaries must contain the key "tag".
    :raises KeyError: Dictionaries must contain the key "r0".
    :raises TypeError: Key "tags" must contain a binary char.
    :raises ValueError: Key "tags" must contain a binary char.
    :raises ValueError: Key "r0" must contain be a list.
    :raises TypeError: Key "r0" must contain be a list.
    :raises ValueError: Key "r0" must contain a list of length two.
    :raises ValueError: Key "r0" must contain a list of two real numbers.
    '''
    if not isinstance(unit_cell, list):
        raise TypeError('\n\nParameter unit_cell must be a list.\n')
    for dic in unit_cell:
        if 'tag' not in dic:
            raise KeyError('\n\nDictionaries must contain the key "tag".\n')
        if 'r0' not in dic:
            raise KeyError('\n\nDictionaries must contain the key "r0".\n')
        if not isinstance(dic['tag'], bytes):
            raise TypeError('\n\Key "tags" must contain a binary char.\n')    
        if not len(dic['tag']) == 1:
            raise ValueError('\n\Key "tags" must contain a binary char.\n')    
        if not isinstance(dic['r0'], list):
            raise TypeError('\n\Key "r0" must contain be a list.\n')
        if not len(dic['r0']) == 2:
            raise ValueError('\n\Key "r0" must contain a list of length two.\n')    
        if not isinstance(dic['r0'][0], (int, float)) or not isinstance(dic['r0'][1], (int, float)):
            raise ValueError('\n\Key "r0" must contain a list of two real numbers.\n')

    
def test_prim_vec(prim_vec):
    '''
    Check parameter *prim_vec*.

    :raises TypeError: Parameter prim_vec must be a dictionary.
    :raises KeyError: Parameter prim_vec must contain the key "norm".
    :raises KeyError: Parameter prim_vec must contain the key "angle".
    :raises TypeError: Parameter  prim_vec["angle"] must be a real number.
    :raises TypeError: Parameter prim_vec["norm"] must be a real number.
    :raises ValueError: Parameter prim_vec["norm"] must be a positive number.
    '''
    if not isinstance(prim_vec, dict):
        raise TypeError('\n\nParameter prim_vec must be a dictionary.\n')
    if 'norm' not in prim_vec:
        raise KeyError('\n\nParameter prim_vec must contain the key "norm".\n')
    if 'angle' not in prim_vec:
        raise KeyError('\n\nParameter prim_vec must contain the key "angle".\n')
    if not isinstance(prim_vec['angle'], (int, float)):
        raise TypeError('\n\nParameter  prim_vec["angle"] must be a real number.\n')
    if not isinstance(prim_vec['norm'], (int, float)):
        raise TypeError('\n\nParameter prim_vec["norm"] must be a real number.\n')
    if prim_vec['norm'] <= 0:
        raise ValueError('\n\nParameter prim_vec["norm"] must be a positive number.\n')
    if prim_vec['norm'] <= 0.2:
        raise ValueError('\n\nParameter prim_vec["norm"] must be a larger than 0.2.\n')


def test_get_lattice(n1, n2):
    '''
    Check method *lattice*.

    :raises TypeError: Parameter n1 must be an integer.
    :raises TypeError: Parameter n2 must be an integer.
    :raises ValueError: Parameter n1 must be a positive integer.
    :raises ValueError: Parameter n2 must be a positive integer.
    '''
    if not isinstance(n1, int):
        raise TypeError('\n\nParameter n1 must be an integer.\n')
    if not isinstance(n2, int):
        raise TypeError('\n\nParameter n2 must be an integer.\n')
    if n1 < 1:
        raise ValueError('\n\nParameter n1 must be a positive integer.\n')
    if n2 < 1:
        raise ValueError('\n\nParameter n2 must be a positive integer.\n')


def test_coor(coor):
    '''
    Check if *get_lattice* has been called (*self.coor* not empty).
    :raises RuntimeError: Run method get_lattice first.
    '''
    if coor.size == 0:
        raise RuntimeError('\n\nRun method get_lattice first.\n')


def test_remove_sites(index, sites):
    '''
    Check method *remove_sites*.

    :raises TypeError: Parameter index must be a list.
    :raises ValueError: Parameter index must be a list of integers.
    :raises ValueError: Indices must be between 0 and sites -1.
      of integers between 0 and sites
    '''
    if not isinstance(index, list):
        raise TypeError('\n\nParameter index must be a list.\n')
    if not all(isinstance(i, int) for i in index):
        raise ValueError('\n\nParameter index must be a list of integers.\n')
    if not all(-1 < i < sites for i in index):
        raise ValueError('\n\nElements of index must be between 0 and sites - 1.\n')


def test_limits(rlim):
    '''
    Check methods *remove_sites_x* and *remove_sites_y*.

    :raises TypeError: Parameter n1 must be an integer.
    :raises TypeError: Parameter n2 must be an integer.
    :raises ValueError: Parameter n1 must be a positive integer.
    :raises ValueError: Parameter n2 must be a positive integer.
    :raises ValueError: Parameter n2 must be a positive integer.
    '''
    if not isinstance(rlim, list):
        raise TypeError('\n\nParameter xlim/ylim must be a list.\n')
    if not len(rlim) == 2:
        raise TypeError('\n\nParameter xlim/ylim must be a list of length two.\n')
    if not isinstance(rlim[0], (int, float)) or not isinstance(rlim[1], (int, float)):
        raise TypeError('\n\nParameter xlim/ylim must be composed of two real numbers.\n')
    if not rlim[0] < rlim[1]:
        raise ValueError('\n\nFirst element of xlim/ylim must be smaller than the second element.\n')


"""
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
"""


class latticeTB(object):
    '''
    Build up 1D or 2D lattice.
    In 2D, a lattice can be expressed by: 
    :math:`\mathbf{R} = n_1\mathbf{a}_1 + n_2\mathbf{a}_2`
    where :math:`\mathbf{a}_1` and `\mathbf{a}_2` are the two primitive vectors.
    '''
    def __init__(self, unit_cell, prim_vec):
        '''
        :param unit_cell: List of dictionaries. One dictionary per site within the unit cell.
          Each dictionary with two keys: 
            * 'tag', label of the site.
            * 'r0', position. 
        :param prim_vec: Dictionary. Define the two primitive vectors.
           Dictionary with two keys: 
            * 'norm', norm of both primitive vectors: `|\mathbf{a}_1|=|\mathbf{a}_2| = norm` 
            * 'angle', angle between the two primitive:
                * :math:`\mathbf{a}_1| = norm (1, 0)` .
                * :math:`|\mathbf{a}_1|= norm (\cos angle, \sin \angle )`.
        '''
        test_unit_cell(unit_cell)
        test_prim_vec(prim_vec)
        self.unit_cell = unit_cell
        self.prim_vec = prim_vec
        self.n1, self.n2 = None, None  # number of unit cells along a1 and a2
        self.coor = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        self.sites = 0  # Site number

    def get_lattice(self, n1, n2=1):
        '''
        Get the lattice positions

        :param n1: Number of unit cells along :math:`x`.
        :param n2: Default value 1. Number of unit cells along :math:`y`.
        '''
        test_get_lattice(n1, n2)
        self.coor = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        self.n1, self.n2 = n1, n2
        x = self.prim_vec['norm'] * np.arange(n1, dtype='f16')
        y = np.sin(PI / 180. * self.prim_vec['angle']) * self.prim_vec['norm'] * np.arange(n2, dtype='f16')
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
        self.coor['x'] = self.coor['x'].round(6)
        self.coor['y'] = self.coor['y'].round(6)
        self.sites = len(self.coor['tag'])

    def remove_sites(self, index):
        '''
        Remove sites defined by their indices
        (use method lattice(plt_index=True) of plotTB 
        to get access to the site indices).

        :param index: List of site indices to be removed.
        '''
        test_coor(self.coor)
        test_remove_sites(index, self.sites)
        mask = np.ones(self.sites, bool)
        mask[index] = False
        self.coor = self.coor[mask]
        self.sites = self.coor.size

    def remove_dangling(self):
        '''
        Remove dangling sites
        (sites connected with just another site).
        '''
        test_coor(self.coor)
        while True:
            dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
            dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
            dis = np.sqrt(dif_x ** 2 + dif_y ** 2).round(2)
            dis_unique = np.unique(dis)
            len_hop = dis_unique[1]
            ind = np.argwhere((dis > len_hop-1e-3) & (dis < len_hop+1e-3))
            dang = []
            for i in range(self.sites):
                if (ind[:, 0] == i).sum() == 1:
                    dang.append(i)
            self.coor = np.delete(self.coor, dang, axis=0)
            self.sites -= len(dang)
            if dang == []:
                break

    def remove_sites_x(self, xlim):
        '''
        Remove site not in xlim.

        :param xlim: List. xlim = [xmin, xmax]. 
        '''
        test_coor(self.coor)
        test_limits(xlim)
        ind = (self.coor['x'] >= xlim[0]-0.01) & (self.coor['x'] < xlim[1]+0.01)
        ind = (self.coor['x'] >= xlim[0]-0.01) & (self.coor['x'] < xlim[1]+0.01)
        self.coor = self.coor[ind]
        self.sites = len(self.coor)

    def remove_sites_y(self, ylim):
        '''
        Remove sites not in ylim.

        :param xlim: List. xlim = [xmin, xmax]. 
        '''
        test_coor(self.coor)
        test_limits(ylim)
        ind = (self.coor['y'] >= ylim[0]-0.01) & (self.coor['y'] < ylim[1]+0.01)
        self.coor = self.coor[ind]
        self.sites = len(self.coor)

    def coor_center(self):
        '''
        Fix the center of mass of the lattice at (0, 0).
        '''
        test_coor(self.coor)
        self.coor['x'] -= np.mean(self.coor['x'])
        self.coor['y'] -= np.mean(self.coor['y'])
