'''
    The **chainTB** module models a one-dimensional tight-binding dimer chain
    with real and/or complex valued hoppings and onsite energies.
    
    **chainTB** is available at https://github.com/cpoli/TB
'''

from latticeTB import *
from plotTB import *
from eigTB import *
from propTB import *
from math import pi


class eigChain(eigTB):
    '''
    Child of the class **eigTB**. Dedicated to chains along :math:`x`. 
    '''
    def __init__(self, lat):
        eigTB.__init__(self, lat)

    def set_hop_alt(self, ho):
        '''
        Set chain alternating hoppings.

        :param ta: Hopping :math:`t_{ab}`.
        :param tb: Default value 0. Hopping :math:`t_{ba}`.  
            If tb=0, t_{ba}=t_{ab}
        '''
        bonds = self.sites - 1
        pos_x = self.coor['x']
        dif_x = pos_x.reshape(self.sites, 1) - pos_x
        ind = np.argwhere(dif_x == - 1)
        hop = np.zeros(bonds, dtype={'names':['i', 'j', 't', 'tag'],
                                                         'formats':['i4','i4','c16', 'S2']})
        hop['i'][:] = ind[:, 0]
        hop['j'][:] = ind[:, 1]
        hop['t'][0::2] = ho[0]
        hop['t'][1::2] = ho[1]
        hop['tag'][0::2] = 'ab'
        hop['tag'][1::2] = 'ba'
        self.hop = hop

    def set_dim_defect(self, dim):
        '''
        Set dimerization defects.

        :param dim: Array. Indices of the dimerization defects.
        '''
        if not self.hop['t'].any():
            raise Exception('\n\nRun method get_hop() first.\n')
        dim = np.array([dim])
        for d in dim:
            t1 = self.hop['t'][d]
            t2 = self.hop['t'][d-1]
            self.hop['t'][d::2] = t2
            self.hop['t'][d+1::2] = t1
            t1 = self.hop['t'][d-1]
            t2 = self.hop['t'][d]

    def get_params(self):
        '''
        Set parameters used to store figures.

        :param on:  Array. Sublattices onsite energies.
        '''
        return {'d_h': self.alpha, 'd_o': self.alpha_on}

class plotChain(plotTB):
    '''
    Child of class **plotTB**. Dedicated to Graphene. 
    '''
    def __init__(self, sys):
            plotTB.__init__(self, sys)

    def plt_dispersion(self, ta, tb=0., lw=5, fs=20):
        '''
        Plot  infinite chain dispersion relation.

        :param ta: Hopping :math:`t_{ab}`.
        :param tb: Default value 0. Hopping :math:`t_{ba}`.    
        :param lw: Default value 5. Linewidth.
        :param fs: Default value 20. Fontsize.
        '''
        if tb == 0:
            tb = ta
        k = np.linspace(-pi, pi, 100)
        e = np.sqrt(ta **2 + tb **2 + 2 * ta * tb * np.cos(k))
        fig, ax = plt.subplots()
        plt.xlim([-pi, pi])
        plt.ylim([-(ta+tb+.2), ta+tb+.2])
        ax.set_xticklabels(['$\pi$', '$\pi/2$', '$0$', '$\pi/2$', '$\pi$'], fontsize=2*fs/3)
        ax.set_xticks([-pi, -pi/2, 0, pi/2, pi])
        plt.plot(k, e, 'b', lw=lw, label='coducting band')
        plt.plot(k, -e, 'r', lw=lw, label='valence band')
        plt.xlabel('$k$', fontsize=fs)
        plt.ylabel('$E(k)$', fontsize=fs)
        plt.title('$E(k)$', fontsize=fs)
        plt.legend(loc='center')
        plt.draw()
        return fig

    def plt_chain_hop(self, ms=20, fs=20):
        '''
        Plot  chain in hopping space.

        :param ms: Default value 20. Markersize.
        :param fs: Default value 20. Font size.
        '''
        coor_hop = np.zeros(self.sys.sites, {'names':['x', 'y', 'tag'], 'formats':['f8', 'f8', 'S1']})
        coor_hop['x'][1:] = np.cumsum(self.sys.hop['t'][:self.sys.sites-1]).real
        coor_hop['tag'] = self.sys.coor['tag']
        fig, ax = plt.subplots(figsize=(10, 2.8))
        plt.title('Lattice in hopping space', fontsize=fs)
        for t, c in zip(self.sys.tags, self.colors):
                plt.plot(coor_hop['x'][coor_hop['tag'] == t],
                            coor_hop['y'][coor_hop['tag'] == t],
                            'o', color=c, ms=ms, markeredgecolor='none')
        plt.xlim([np.min(coor_hop['x'])-1.5, np.max(coor_hop['x'])+1.5])
        plt.ylim([np.min(coor_hop['y'])-1.5, np.max(coor_hop['y'])+1.5])
        plt.xticks([])
        plt.yticks([])
        plt.draw()
        return fig
