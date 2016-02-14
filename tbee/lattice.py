import numpy as np
import matplotlib.pyplot as plt
import tbee.error_handling as error_handling


PI = np.pi


#################################
# CLASS LATTICE
#################################


class lattice():
    r'''
    Build up 1D or 2D lattice.
    Lattice is defined by the discrete operation:

    .. math::

        \mathbf{R} = n_1\mathbf{a}_1 + n_2\mathbf{a}_2

    where :math:`\mathbf{a}_1` and :math:`\mathbf{a}_2` are the two primitive
    vectors and :math:`n_1` and :math:`n_2` are the number of unit cells along 
    :math:`\mathbf{a}_1` and :math:`\mathbf{a}_2`.

    :param unit_cell: List of dictionaries. 
     One dictionary per site within the unit cell. Each dictionary has two keys:

        * 'tag', Binary Char. label of the associated sublattice.
        * 'r0', Tuple. Position.
    :param prim_vec: List of tuples. 
     Define the primitive vectors. List of one/two tuples for 1D/2D respectively:
     
        * Tuple, cartesian coordinate of the primitive vector :math:`\mathbf{a}_1`.
        * Tuple, cartesian coordinate of the primitive vector :math:`\mathbf{a}_2`.

    Example usage::

        # Line-Centered Square lattice
        unit_cell = [{'tag': b'a', r0=(0., 0.)}, {'tag': b'a', r0=(0., 1.)}]
        prim_vec = [(0, 2), (2, 0)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    '''

    def __init__(self, unit_cell, prim_vec):
        error_handling.unit_cell(unit_cell)
        error_handling.prim_vec(prim_vec)
        self.unit_cell = unit_cell
        self.prim_vec = prim_vec
        self.tags = np.unique(np.array([dic['tag'] for dic in self.unit_cell]))
        self.n1, self.n2 = 0, 0
        self.coor = np.array([], dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'S1')])
        self.sites = 0

    def get_lattice(self, n1, n2=1):
        '''
        Get the lattice positions.

        :param n1: Positive Integer.
            Number of unit cells along :math:`\mathbf{a}_1`.
        :param n2: Positive Integer. Default value 1.
            Number of unit cells along :math:`\mathbf{a}_2`.

        Example usage::

            # Line-Centered Square lattice
            unit_cell = [{'tag': b'a', r0=(0., 0.)}, {'tag': b'a', r0=(0., 1.)}]
            prim_vec = [(0, 2), (2, 0)]
            lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
            lat.get_lattice(n1=4, n2=5)
        '''
        error_handling.get_lattice(self.prim_vec, n1, n2)
        sites_uc = len(self.unit_cell)
        sites_tag = n1*n2
        self.sites = sites_uc * sites_tag
        self.coor = np.empty(self.sites, dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'S1')])
        self.n1, self.n2 = n1, n2
        x = self.prim_vec[0][0] * np.arange(n1, dtype='f8')
        y = self.prim_vec[0][1] * np.arange(n1, dtype='f8')
        xx = np.empty(n1*n2)
        yy = np.empty(n1*n2)
        xx[:n1] = x
        yy[:n1] = y
        for i in range(1, n2):
            xx[i*n1: (i+1)*n1] = x + i * self.prim_vec[1][0]
            yy[i*n1: (i+1)*n1] = y + i * self.prim_vec[1][1]
        for i, dic in enumerate(self.unit_cell):
            self.coor['x'][i*sites_tag: (i+1)*sites_tag] = xx + dic['r0'][0]
            self.coor['y'][i*sites_tag: (i+1)*sites_tag] = yy + dic['r0'][1]
            self.coor['tag'][i*sites_tag: (i+1)*sites_tag] = dic['tag']
        self.coor = np.sort(self.coor, order=('y', 'x'))

    def add_sites(self, coor):
        '''
        Add sites.

        :param coor: Structured array with keys: {'x', 'y', 'tag'}.

        Example usage::

            # Square lattice
            unit_cell = [{'tag': b'a', r0=(0., 0.)}]
            prim_vec = [(0, 1), (1, 0)]
            lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
            lat.get_lattice(n1=2, n2=2)
            coor = np.array([(-1., -1, b'b'), (-2., -2, b'c')], 
                                      dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'S1')])
            lat.add_sites(coor)
        '''
        error_handling.coor(coor)
        self.coor = np.concatenate([self.coor, coor])
        self.sites += len(coor)
        self.tags = np.unique(np.concatenate([self.tags, coor['tag']]))
        self.coor = np.sort(self.coor, order=('y', 'x'))

    def remove_sites(self, index):
        '''
        Remove sites defined by their indices
        (use method lattice.plot(plt_index=True)
        to get access to the site indices).

        :param index: List. Site indices to be removed.

        Example usage::

            # Square lattice
            unit_cell = [{'tag': b'a', r0=(0., 0.)}]
            prim_vec = [(0, 1), (1, 0)]
            lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
            lat.get_lattice(n1=2, n2=2)
            lat.remove_sites([0, 2])
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

    def ellipse_in(self, rx, ry, x0, y0):
        '''
        Select sites according to 

        .. math:: 

            (x-x_0)^2/a^2+(y-y_0)^2/b^2 < 1\,  .

        :param list_hop: List of Dictionary (see set_hopping definition).
        :param rx: Positive Real number. Radius along :math:`x`. 
        :param ry: Positive Real number. Radius along :math:`y`.
        :param x0: Real number. :math:`x` center. 
        :param y0: Real number. :math:`y` center.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.positive_real(rx, 'rx')
        error_handling.positive_real(ry, 'ry')
        error_handling.real_number(x0, 'x0')
        error_handling.real_number(y0, 'y0')
        self.coor = self.coor[(self.coor['x'] -x0) ** 2 / rx ** 2 +  \
                                        (self.coor['y'] -y0) ** 2 / ry ** 2 < 1.]
        self.sites = len(self.coor)

    def ellipse_out(self, rx, ry, x0, y0):
        '''
        Select sites according to

        .. math:: 

            (x-x_0)^2/a^2+(y-y_0)^2/b^2 > 1\,  .


        :param list_hop: List of Dictionary (see set_hopping definition).
        :param rx: Positive Real number. Radius along :math:`x`. 
        :param ry: Positive Real number. Radius along :math:`y`.
        :param x0: Real number. :math:`x` center. 
        :param y0: Real number. :math:`y` center.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.positive_real(rx, 'rx')
        error_handling.positive_real(ry, 'ry')
        error_handling.real_number(x0, 'x0')
        error_handling.real_number(y0, 'y0')
        self.coor = self.coor[(self.coor['x'] -x0) ** 2 / rx ** 2 +  \
                                        (self.coor['y'] -y0) ** 2 / ry ** 2 > 1.]
        self.sites = len(self.coor)

    def center(self):
        '''
        Fix the center of mass of the lattice at (0, 0).
        '''
        error_handling.empty_coor(self.coor)
        self.coor['x'] -= np.mean(self.coor['x'])
        self.coor['y'] -= np.mean(self.coor['y'])

    def rotation(self, theta):
        r'''
        Rotate the lattice structure by the angle :math:`\theta`.

        :param theta: Rotation angle in degrees.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.real_number(theta, 'theta')
        theta *= PI / 360
        for dic in self.unit_cell:
            x  = self.coor['x'] - dic['r0'][0]
            y  = self.coor['y'] - dic['r0'][1]
            self.coor['x'] = x * np.cos(theta) - y * np.sin(theta) + dic['r0'][0]
            self.coor['y'] = y * np.cos(theta) + x* np.sin(theta) + dic['r0'][1]

    def clean_coor(self):
        '''
        Keep only the sites with different coordinates.
        '''
        error_handling.empty_coor(self.coor)
        coor = self.coor[['x', 'y']].copy()
        coor['x'], coor['y'] = self.coor['x'].round(4), self.coor['y'].round(4)
        _, idx = np.unique(coor, return_index=True)
        self.coor = self.coor[idx]
        self.sites = len(self.coor)

    def __add__(self, other):
        '''
        Overloading operator +.
        '''
        error_handling.lat(other)
        error_handling.empty_coor(self.coor)
        error_handling.empty_coor(other.coor)
        coor = np.concatenate([self.coor, other.coor])
        tags = np.concatenate([self.tags, other.tags])
        lat = lattice(unit_cell=self.unit_cell, prim_vec=self.prim_vec)
        lat.add_sites(coor)
        lat.sites = self.sites + other.sites
        lat.tags = np.unique(tags)
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

        .. note::

            The tags are not considered in the lattice subtraction.
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

        .. note::

            The tags are not considered in the lattice subtraction.
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

    def plot(self, ms=20, fs=20, plt_index=False, axis=False, figsize=None):
        '''
        Plot lattice in hopping space.

        :param ms: Positive number. Default value 20. Markersize.
        :param fs: Positve number. Default value 20. Fontsize.
        :param plt_index: Boolean. Default value False. Plot site labels.
        :param axis: Boolean. Default value False. Plot axis.
        :param figsize: Tuple. Default value None. Figsize.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_coor(self.coor)
        error_handling.positive_real(ms, 'ms')
        error_handling.positive_real(fs, 'fs')
        error_handling.boolean(plt_index, 'plt_index')
        error_handling.boolean(axis, 'axis')
        if figsize is None:
            figsize = (5, 5)
        error_handling.list_tuple_2elem(figsize, 'figsize')
        error_handling.positive_real(figsize[0], 'figsize[0]')
        error_handling.positive_real(figsize[1], 'figsize[1]')
        fig, ax = plt.subplots(figsize=figsize)
        # plot sites
        colors = ['b', 'r', 'g', 'y', 'm', 'k']
        for color, tag in zip(colors, self.tags):
            plt.plot(self.coor['x'][self.coor['tag'] == tag],
                        self.coor['y'][self.coor['tag'] == tag],
                       'o', color=color, ms=ms, markeredgecolor='none')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(self.coor['x'])-1., np.max(self.coor['x'])+1.])
        ax.set_ylim([np.min(self.coor['y'])-1., np.max(self.coor['y'])+1.])
        if not axis:
            ax.axis('off')
        # plot indices
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sites)]
            for l, x, y in zip(indices, self.coor['x'], self.coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                                    textcoords='offset points',
                                    ha='right', va='bottom', size=fs)
        plt.draw()
        return fig

    def show(self):
        """
        Emulate Matplotlib method plt.show().
        """
        plt.show()