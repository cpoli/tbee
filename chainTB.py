'''
    The **chainTB** package models a one-dimensional tight-binding dimer chain
    with real and/or complex valued hoppings and onsite energies.
    Defects can be introduced by locally modifying the onsite energies
    and/or the hoppings.
    
    **chainTB** is available at https://github.com/cpoli/chainTB
'''

from latticeTB import *
from plotTB import *
from eigTB import *
from propagationTB import *
from math import pi
import numpy as np


class eigChain(eigTB):
    '''
    Child of the class **eigTB**. Dedicated to chains along :math:`x`. 
    '''
    def __init__(self, lat):
        eigTB.__init__(self, lat)
        self.c = 0. 

    def set_hop(self, ta, tb=0):
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
        if tb == 0.:
            tb = ta
        hop['i'][:] = ind[:, 0]
        hop['j'][:] = ind[:, 1]
        hop['t'][0::2] = ta
        hop['t'][1::2] = tb
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
        return {'d_h': self.alpha_hop, 'd_o': self.alpha_on}

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
        coor_hop['x'][1:] = np.cumsum(self.sys.hop['t']).real
        coor_hop['tag'] = self.sys.coor['tag']
        fig, ax = plt.subplots(figsize=(10, 2.8))
        plt.title('Lattice in hopping space', fontsize=fs)
        for t, c in zip(self.sys.tags, self.colors):
                plt.plot(coor_hop['x'][coor_hop['tag'] == t],
                            coor_hop['y'][coor_hop['tag'] == t],
                            'o', color=c, ms=ms, markeredgecolor='none')
        xlim = [np.min(coor_hop['x'])-1.5, np.max(coor_hop['x'])+1.5]
        ylim = [np.min(coor_hop['y'])-1.5, np.max(coor_hop['y'])+1.5]
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xticks([])
        plt.yticks([])
        plt.draw()
        return fig

    """
    TO BE IMPLEMENTED
    def plt_geometry(self, hops, dz=.05, steps_pump=100, fs=20):
        '''
        Plot hopping geometry.

        :param hops: List of hopping configuration.
        :param dz: Default value 0.05. Step.
        :param steps_pump: Default value 100. Number of time steps stent in
            one hopping configuration, and to go from one hopping configuration
            to the next one.
        :param fs: Default value 20. Fontsize.
        :param save: Default value False. Save the figure in the directory given 
            by the method *dir_name* with name *geometry*.
        '''
        hop_new = np.array([hop_new])
        if not all(len(hop) == self.bonds for hop in hop_new):
            raise Exception('\n\nArgument hop_new is not composed '
                            'of size (sites-1, N).\n')
        hop_save = self.hop
        num_pump = len(hop_new)
        steps_tot = (num_pump+2)*steps_pump
        geo = np.zeros((self.sites, steps_tot))
        z = dz * np.arange(0, steps_tot)
        # Before pumping
        hop_cum = np.cumsum(self.hop.real)
        pos = np.zeros(self.sites)
        pos[1:] = hop_cum
        x, y = np.meshgrid(np.arange(0, steps_pump), pos)
        geo[:, :steps_pump] = y
        steps = range(0, steps_pump)
        for j in steps:
            geo[1:self.sites, j] = np.cumsum(self.hop.real)
        # Pumping
        trans = np.linspace(0., 1., steps_pump)
        for n in range(num_pump):
            hop_before = self.hop
            steps = range((n+1)*steps_pump, (n+2)*steps_pump)
            for ind, val in enumerate(steps):
                self.hop = \
                    trans[ind]*(hop_new[:][n] - hop_before) + hop_before
                geo[1:self.sites, val] = np.cumsum(self.hop.real)
        # After pumping
        steps = range((1+num_pump)*steps_pump, steps_tot)
        for j in steps:
            geo[1:self.sites, j] = np.cumsum(self.hop.real)
        self.hop = hop_save
        fig, ax = plt.subplots()
        plt.title('HOPPINGS TIME EVOLUTION')
        ax.set_ylim([-1, geo[self.sites-1, 0]+1])
        for i in range(0, self.sites, 2):
            plt.plot(z, geo[i, :], 'b', linewidth=4)
        for i in range(1, self.sites, 2):
            plt.plot(z, geo[i, :], 'r', linewidth=4)
        plt.xlabel('$z$', fontsize=fs)
        plt.xlim(z[0], z[-1])
        frame = plt.gca()
        frame.axes.get_yaxis().set_visible(False)
        return fig
        """