import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, sin, cos, pi, tan, atan2
import sys
from plotTB import *

np.set_printoptions(precision=2)


class latticeTB():
    '''
    Build up a lattice.
    '''
    def __init__(self, tags, ri, nor, ang):
        '''
        :param tags: The tags (binary chars) associated to each sublattice.
        :param ri: Initial positions of sites belonging to the unit cell.
        :param nor: Primitive vectors norm.
        :param ang: Angle between primitive vectors.
        '''
        self.tags = np.array(tags)
        self.nx, self.ny = 0, 0  # no unit cells along x and y
        self.ri = ri  # initial positions of the sites
        self.nor = nor  # primitive vectors norm 
        self.ang = ang  # angle between primitive vectors     
        self.coor = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.sites = 0  # Site number   

    def get_lattice(self, nx, ny):
        '''
        Get the lattice positions
        
        :param nx: Number of unit cells along :math:`x`.
        :param ny: Number of unit cells along :math:`y`.
        '''
        self.coor = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.nx, self.ny = nx, ny
        x = self.nor * np.arange(nx, dtype='f8')
        y = sin(self.ang) * self.nor * np.arange(ny, dtype='f8')
        xx, yy = np.meshgrid(x, y)
        if self.ang:
            xx += yy / tan(self.ang) 
        xx = np.ravel(xx)
        yy = np.ravel(yy)
        coor_tag = np.zeros(len(xx), {'names':['x', 'y', 'tag'], 'formats':['f8', 'f8', 'S1']})
        for i, t in enumerate(self.tags):
            coor_tag['x'] = xx + self.ri[i][0]
            coor_tag['y'] = yy + self.ri[i][1]
            coor_tag['tag'] = t             
            self.coor = np.concatenate([self.coor, coor_tag])
        self.coor =  np.sort(self.coor, order=('y', 'x'))
        if (self.ang == pi/2) and (int(self.nor) == self.nor):
            self.coor['x'] = np.rint(self.coor['x'])
            self.coor['y'] = np.rint(self.coor['y'])
        self.sites = len(self.coor['tag'])

    def remove_sites(self, ind):
        '''
        Remove sites define their indices
        
        :param ind: Array of site indices to be removed.
        '''
        mask = np.ones(self.sites, bool)
        mask[ind] = False
        self.coor = self.coor[mask] 
        self.sites = len(self.coor['tag'])

    def remove_dangling(self, nor_bond):
        '''
        Remove dangling sites
        (sites connected with just another site).
        
        :param nor_bond: Norm of the bonds to be checked.
        '''
        while True: 
            dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
            dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
            dis = np.sqrt(dif_x **2 + dif_y **2)
            ind = np.argwhere((dis > nor_bond-.1) & (dis < nor_bond+.1))
            dang =[]
            for i in range(self.sites):
                if (ind[:, 0] == i).sum() == 1:
                    dang.append(i)
            self.coor = np.delete(self.coor, dang, axis=0)
            self.sites -= len(dang)
            if dang == []:
                break

    def plt_lattice(self, colors=[], ms=30, fs=20, plt_index=False, figsize=None):
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
        if colors == []:
            colors = ['b', 'r', 'g', 'y', 'm', 'k']
        fig, ax = plt.subplots(figsize=figsize)
        # sites
        for c, t in zip(colors, self.tags):
            str_tag = t.decode('ascii').upper()
            plt.plot(self.coor['x'][self.coor['tag'] == t],
                        self.coor['y'][self.coor['tag'] == t], 'o', color=c,
                        ms=ms, label=str_tag, markeredgecolor='none')
        ax.axis('off')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(self.coor['x'])-0.5, np.max(self.coor['x'])+0.5])
        ax.set_ylim([np.min(self.coor['y'])-0.5, np.max(self.coor['y'])+0.5])
        # labels
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sites)]
            for l, x, y in zip(indices, self.coor['x'], self.coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        return fig
