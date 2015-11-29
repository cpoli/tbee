from latticeTB import *
from plotTB import *
from eigTB import *
from propTB import *
from scipy.optimize import fsolve

def test_set_nearest_neighbor_hop(t_ab, t_ba, t_ac, t_ca):
    '''
    Check method *set_nearest_neighbor_hop*.

    :raises TypeError: Parameter t_ab must be a number.
    :raises TypeError: Parameter t_ba must be a number.
    :raises TypeError: Parameter t_ac must be a number.
    :raises TypeError: Parameter t_ca must be a number.
    '''
    if not isinstance(t_ab, (int, float, complex)):
        raise TypeError('\n\nParameter t_ab must be a number.\n')
    if not isinstance(t_ba, (int, float, complex)):
        raise TypeError('\n\nParameter t_ba must be a number.\n')
    if not isinstance(t_ac, (int, float, complex)):
        raise TypeError('\n\nParameter t_ac must be a number.\n')
    if not isinstance(t_ca, (int, float, complex)):
        raise TypeError('\n\nParameter t_ca must be a number.\n')


def test_set_next_nearest_neighbor_hop(t_ab, t_ba, t_ac, t_ca, c):
    '''
    Check method *set_nearest_neighbor_hop*.

    :raises TypeError: Parameter t_ab must be a number.
    :raises TypeError: Parameter t_ba must be a number.
    :raises TypeError: Parameter t_ac must be a number.
    :raises TypeError: Parameter t_ca must be a number.
    :raises TypeError: Parameter c must be a real number.
    '''
    test_set_nearest_neighbor_hop(t_ab, t_ba, t_ac, t_ca)
    if not isinstance(c, (int, float)):
        raise TypeError('\n\nParameter t_ca must be a number.\n')

def test_set_defect_dimer(r):
    '''
    Check method *set_defect_dimer_x* and *set_set_defect_dimer_y*.

    :raises TypeError: Parameter x/y must be a real number.
    :raises ValueError: Parameter x/y must be a positive number.
    '''
    if not isinstance(r, (int, float)):
        raise TypeError('\n\nParameter x (or y) must be a real number.\n')
    if r < 0:
        raise ValueError('\n\nParameter x (or y) must be a positive number.\n')


def test_set_disorder(alpha):
    '''
    Check methods *set_disorder*.

    :raises TypeError: Parameter alpha must be a real number.
    :raises ValueError: Parameter n must be a positive number.
    '''
    if not isinstance(alpha, (int, float)):
        raise TypeError('\n\nParameter alpha must be a real number.\n')
    if alpha < 0:
        raise ValueError('\n\nParameter alpha must be a positive number.\n')


class liebTB(latticeTB):
    '''
    Child of the class **latticeTB**. Dedicated to the Lieb lattice. 
    '''
    def __init__(self):
        ri = [[0., 0.], [1., 0.], [0., 1.]]
        tags = [b'a', b'b', b'c']
        nor = 2.
        ang = 90.
        latticeTB.__init__(self, ri=ri, tags=tags, nor=nor, ang=ang)


class liebEig(eigTB):
    def __init__(self, lat):
        eigTB.__init__(self, lat)
        # nearest neighbor hopping
        self.t_ab, self.t_ba, self.t_ac, self.t_ca = None, None, None, None
        # next nearest neighbor hopping
        self.t_pp, self.t_pm, self.t_mp, self.t_mm = None, None, None, None
        self.c = None
        # dimerization defect bool
        self.defect_dimer_x, self.defect_dimer_y = None, None
        self.coor_hop = np.array([], dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])

    def set_nearest_neighbor_hop(self, t_ab, t_ba, t_ac, t_ca):
        test_set_nearest_neighbor_hop(t_ab, t_ba, t_ac, t_ca)
        print('ab', t_ab, 'ba', t_ba, 'ac', t_ac, 'ca', t_ca)
        self.set_hop([{'n': 1, 'hop': [{'tag': b'ab', 't': t_ab}, {'tag': b'ba', 't': t_ba},
                                               {'tag': b'ac', 't': t_ac}, {'tag': b'ca', 't': t_ca}]}])
        self.t_ab, self.t_ba, self.t_ac, self.t_ca = t_ab, t_ba, t_ac, t_ca

    def set_next_nearest_neighbor_hop(self, c):
        test_set_nearest_neighbor_hop(self.t_ab, self.t_ba, self.t_ac, self.t_ca)
        self.set_hop([{'n': 2, 'hop': [{'tag': b'bc', 't': 0}, {'tag': b'cb', 't': 0}]}])
        self.rename_hop_tag([{'tag_new': b'++', 'n': 2, 'ang': 45, 'tag': b'bc'},
                                       {'tag_new': b'-+', 'n': 2, 'ang': 135, 'tag': b'bc'}, 
                                       {'tag_new': b'+-', 'n': 2, 'ang': 135, 'tag': b'cb'}, 
                                       {'tag_new': b'--', 'n': 2, 'ang': 45, 'tag': b'cb'}])
        t_pp = c * np.sqrt(self.t_ba + self.t_ac)
        t_pm = c * np.sqrt(self.t_ba + self.t_ca)
        t_mp = c * np.sqrt(self.t_ab + self.t_ac)
        t_mm = c * np.sqrt(self.t_ab + self.t_ca)
        print('next_nearest hoppings')
        print('t++', t_pp, 't--', t_mm, 't+-', t_pm, 't-+', t_mp)
        self.set_hop_with_tag([{'tag': b'++', 't': t_pp}, {'tag': b'+-', 't': t_pm}, 
                                        {'tag': b'-+', 't': t_mp}, {'tag': b'--', 't': t_mm}])
        self.t_pp, self.t_pm, self.t_mp, self.t_mm = t_pp, t_pm, t_mp, t_mm
        self.c = c

    def set_defect_dimer_x(self, x):
        test_set_defect_dimer(x)
        if self.c:
            list_hop=[{'tag': b'ab', 't': self.t_ba}, {'tag': b'ba', 't': self.t_ab},
                          {'tag': b'++', 't': self.t_mp}, {'tag': b'-+', 't': self.t_pp},
                          {'tag': b'--', 't': self.t_pm}, {'tag': b'+-', 't': self.t_mm}]
        else:
            list_hop=[{'tag': b'ab', 't': self.t_ba}, {'tag': b'ba', 't': self.t_ab}]
        self.set_defect_dimer(list_hop=list_hop, x_bottom_left=x, y_bottom_left=0)
        if self.defect_dimer_x:
            raise RuntimeError('\n\n*set_defect_dimer_x* already called')
        self.defect_dimer_x = x
        if isinstance(self.defect_dimer_y, int) and self.c:
            list_hop = [{'tag': b'++', 't': self.t_mm}, {'tag': b'-+', 't': self.t_pm},
                            {'tag': b'--', 't': self.t_pp}, {'tag': b'+-', 't': self.t_mp}]
            self.set_defect_dimer(list_hop=list_hop,
                                                   x_bottom_left=self.defect_dimer_x,
                                                   y_bottom_left=self.defect_dimer_y)

    def set_defect_dimer_y(self, y):
        test_set_defect_dimer(y)
        if self.c:
            list_hop=[{'tag': b'ac', 't': self.t_ca}, {'tag': b'ca', 't': self.t_ac},
                          {'tag': b'++', 't': self.t_pm}, {'tag': b'-+', 't': self.t_mm},
                          {'tag': b'--', 't': self.t_mp}, {'tag': b'+-', 't': self.t_pp}]
        else:
            list_hop=[{'tag': b'ac', 't': self.t_ca}, {'tag': b'ca', 't': self.t_ac}]
        self.set_defect_dimer(list_hop=list_hop, x_bottom_left=0, y_bottom_left=y)
        if self.defect_dimer_y:
            raise RuntimeError('\n\n*set_defect_dimer_y* already called')
        self.defect_dimer_y = y
        if isinstance(self.defect_dimer_x, int) and self.c:
            list_hop = [{'tag': b'++', 't': self.t_mm}, {'tag': b'-+', 't': self.t_pm},
                            {'tag': b'--', 't': self.t_pp}, {'tag': b'+-', 't': self.t_mp}]
            self.set_defect_dimer(list_hop=list_hop,
                                          x_bottom_left=self.defect_dimer_x,
                                          y_bottom_left=self.defect_dimer_y)

    def set_disorder_uniform(self, alpha):
        '''
        Set a non generic disorder.
        Disorder uniform along math:`y` for math:`t_{ab}`  and math:`t_{ba}`,  
        and uniform along math:`x` for math:`t_{ac}`  and math:`t_{ca}`.

        :param alpha: Stength of the disorder.

        .. note ::

            This disorder preserves the zero mode.

        '''
        test_set_disorder(alpha)
        sites_x = len(self.lat.coor['tag'][self.lat.coor['y'] < 1e-4]) 
        sites_y = len(self.lat.coor['tag'][self.lat.coor['x'] < 1e-4]) 
        rand_x = 1. + alpha * rand.uniform(-.5, .5, sites_x-1)
        rand_y = 1. + alpha * rand.uniform(-.5, .5, sites_y-1)
        for y in range(0, sites_y, 2):
            self.hop['t'][((self.hop['tag'] == b'ab') | (self.hop['tag'] == b'ba')) &
                             (self.lat.coor['y'][self.hop['i']] == y)] *= rand_x
        for x in range(0, sites_x, 2):
            self.hop['t'][((self.hop['tag'] == b'ac') | (self.hop['tag'] == b'ca')) &
                             (self.lat.coor['x'][self.hop['i']] == x)] *= rand_y
        self.alpha = alpha

    def set_disorder_pair(self, alpha):
        '''
        Set a non generic disorder.

        :param alpha: Stength of the disorder.

        .. note ::

            This disorder preserves the zero mode.
        '''
        test_set_disorder(alpha)
        sites_x = len(self.lat.coor['tag'][self.lat.coor['y'] < 1e-4]) 
        sites_y = len(self.lat.coor['tag'][self.lat.coor['x'] < 1e-4]) 
        for y in range(0, sites_y, 2):
            rand_x = 1 + alpha * rand.uniform(-.5, .5, sites_x-1)
            rand_x[1::2] = rand_x[0::2]
            self.hop['t'][((self.hop['tag'] == b'ab') | (self.hop['tag'] == b'ba')) & 
		       (self.lat.coor['y'][self.hop['i']] == y)] *= rand_x
        for x in range(0, sites_x, 2):
            rand_y = 1 + alpha * rand.uniform(-.5, .5, sites_y-1)
            rand_y[1::2] = rand_y[0::2]
            self.hop['t'][((self.hop['tag'] == b'ac') | (self.hop['tag'] == b'ca')) &
		       (self.lat.coor['x'][self.hop['i']] == x)] *= rand_y
        self.alpha = alpha
    
    def set_disorder_placket(self, alpha):
        '''
        Set a non generic disorder.

        :param alpha: Stength of the disorder.

        .. note ::

            This disorder preserves the zero mode.

        '''
        self.set_disorder_hop(alpha)
        list_placket = self.get_placket()
        for placket in list_placket:
            t1 = self.hop['t'][(self.hop['i'] == placket[0]) & (self.hop['j'] == placket[1])]
            t2 = self.hop['t'][(self.hop['i'] == placket[1]) & (self.hop['j'] == placket[2])]
            t3 = self.hop['t'][(self.hop['i'] == placket[2]) & (self.hop['j'] == placket[4])]
            t4 = self.hop['t'][(self.hop['i'] == placket[4]) & (self.hop['j'] == placket[7])] 
            t5 = self.hop['t'][(self.hop['i'] == placket[6]) & (self.hop['j'] == placket[7])] 
            t6 = self.hop['t'][(self.hop['i'] == placket[5]) & (self.hop['j'] == placket[6])] 
            t7 = self.hop['t'][(self.hop['i'] == placket[3]) & (self.hop['j'] == placket[5])] 
            t8 = self.hop['t'][(self.hop['i'] == placket[0]) & (self.hop['j'] == placket[3])] 
            t6 = t1 * t3 * t5 * t7 / (t2 * t4 * t8)
            self.hop['t'][(self.hop['i'] == placket[5]) & (self.hop['j'] == placket[6])] = t6
            
    def get_placket(self):
            lx, ly = 2., 2.
            list_placket = []
            for i in range(self.lat.sites):
                ind_placket = np.where((self.lat.coor['x'] - self.lat.coor['x'][i] > - 0.1) & 
                                                    (self.lat.coor['x'] - self.lat.coor['x'][i] < lx + 0.1) &
                                                    (self.lat.coor['y'] - self.lat.coor['y'][i] > - 0.1) & 
                                                    (self.lat.coor['y'] - self.lat.coor['y'][i] < ly + 0.1))
                ind_placket = np.ravel(ind_placket)
                if len(ind_placket) == 8:
                    list_placket.append(ind_placket)
            return list_placket

    def get_coor_hop_dis(self):
        self.get_coor_hop()
        list_placket = self.get_placket()
        for placket in list_placket:
            t3 = self.hop['t'][(self.hop['i'] == placket[2]) & (self.hop['j'] == placket[4])]
            t4 = self.hop['t'][(self.hop['i'] == placket[4]) & (self.hop['j'] == placket[7])] 
            t5 = self.hop['t'][(self.hop['i'] == placket[6]) & (self.hop['j'] == placket[7])] 
            t6 = self.hop['t'][(self.hop['i'] == placket[5]) & (self.hop['j'] == placket[6])] 
            ta = float((t3 + t4).real)
            tb = float((t5 + t6).real)
            param = (self.coor_hop['x'][placket[5]], self.coor_hop['y'][placket[5]], tb,
                   self.coor_hop['x'][placket[2]], self.coor_hop['y'][placket[2]], ta)
            xb, yb, lb, xa, ya, la = param
            if np.sqrt((xb-xa)**2+(yb-ya)**2) > lb+la:
                print('\n\nWARNING WARNING WARNING')
                print('Lattice construction in hopping space is not possible.')
                raise Exception('\n\nRun again or decrease the strength '
                'of the disorder.\n\n')
            angb, anga  = fsolve(equations, ( 0.5*PI, 0.), (param,))
            self.coor_hop['x'][placket[4]] = self.coor_hop['x'][placket[2]] + \
                                                        float(t3.real) * cos(anga)
            self.coor_hop['y'][placket[4]] = self.coor_hop['y'][placket[2]] + \
                                                        float(t3.real) * sin(anga)
            self.coor_hop['x'][placket[6]] = self.coor_hop['x'][placket[5]] + \
                                                            float(t6.real) * cos(angb)
            self.coor_hop['y'][placket[6]] = self.coor_hop['y'][placket[5]] + \
                                                            float(t6.real) * sin(angb) 
            self.coor_hop['x'][placket[7]] = self.coor_hop['x'][placket[6]] + \
                                                            float(t5.real) * cos(angb) 
            self.coor_hop['y'][placket[7]] = self.coor_hop['y'][placket[6]] + \
                                                             float(t5.real) *sin(angb)

def equations(angle, param):
    '''
    Get the intersection of the points at the extremity of two segments.
    Used to build up the lattice in hopping space.
    '''
    (ang2, ang5) = angle
    x2, y2, l2, x5, y5, l5 = param
    return (x5 + l5 * cos(ang5) - x2 - l2 * cos(ang2),
                y5 + l5 * sin(ang5) - y2 - l2 * sin(ang2))

                
                
if __name__ == '__main__':
    t_ab, t_ba, t_ac, t_ca =  4/3, 2/3, 4/3, 2/3
    t_ab, t_ba, t_ac, t_ca = np.round(4/3, 2), np.round(2/3, 2), np.round(4/3, 2), np.round(2/3, 2)
    c = 0.01
    nx, ny = 5, 5
    s = 100
    alpha = 0.5
    lat = liebTB()
    lat.get_lattice(nx=nx, ny=nx)
    lat.remove_dangling(1)

    sys = liebEig(lat)
    plot = plotTB(sys)
    sys.set_nearest_neighbor_hop(t_ab=t_ab, t_ba=t_ba, t_ac=t_ac, t_ca=t_ca)
    sys.set_next_nearest_neighbor_hop(c=c)
    #sys.get_disorder_coor_hop()
    
    sys.set_defect_dimer_x(4)
    sys.set_defect_dimer_y(4)
    sys.set_disorder_hop(alpha)
    sys.get_ham(complex_transpose=True)
    sys.get_eig(eigenvec=True)

    #sys.set_disorder_pair(alpha)
    #sys.set_disorder_placket(alpha)


    zero_mode = sys.get_state_pola(tag_pola=b'a')
    plot.intensity(zero_mode)
    sys.get_coor_hop()
    #plot.lattice_hop(ms=20)

    sys.hop['t'] = 2*np.exp(-sys.hop['t'])
    if alpha == 0:
        sys.get_coor_hop()
    else:
        sys.get_coor_hop_dis()
    fig_coor_real = plot.lattice_hop()

    from collections import OrderedDict
    save = saveTB(sys=sys, params=OrderedDict([('ta', t_ab), ('tb', t_ba), ('tc', t_ca), ('td', t_ca)]), 
                                         dir_name='lieb_supp', ext='pdf')
    save.fig(fig_coor_real, 'lattice_real_gen2')
    plt.show()