import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=2)
import sys

class latticeTB():
    '''
    Build up a lattice.
    '''
    def __init__(self, tags, ri, pv):
        '''
        :param tags: The tags (binary chars) associated to each sublattice.
        :param ri: Initial positions of sites belonging to the unit cell.
        :param pv: Primitive vectors.
        '''
        self.tags = np.array(tags)
        self.nx, self.ny = 0, 0  # no sites along x and y
        self.ri = ri  # initial positions of the sites[[0, 0], [1, 0], [0, 1]]
        self.pv = pv  # primitive vectors       
        self.coor = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.sites = 0  # Site number   

    def get_lattice(self, nx, ny):
        '''
        Get the lattice positions
        
        :param nx: Number of sites along :math:`x`.
        :param ny: Number of sites along :math:`y`.
        '''
        self.coor = np.array([], {'names': ['x', 'y', 'tag'], 'formats': ['f8', 'f8', 'S1']})
        self.nx, self.ny = nx, ny
        for i, t in enumerate(self.tags):
            if self.pv[i][0] == 0:
                self.pv[i][0] = 1
            if self.pv[i][1] == 0:
                self.pv[i][1] = 1
            x_tag = np.arange(self.ri[i][0], self.nx, self.pv[i][0])
            y_tag = np.arange(self.ri[i][1], self.ny, self.pv[i][1])
            x_tag, y_tag = np.meshgrid(x_tag, y_tag)
            x_tag, y_tag = np.ravel(x_tag), np.ravel(y_tag)
            coor_tag = np.zeros(len(x_tag), {'names':['x', 'y', 'tag'], 'formats':['f8', 'f8', 'S1']})
            coor_tag['x'] = x_tag
            coor_tag['y'] = y_tag
            coor_tag['tag'] = t
            self.coor = np.concatenate([self.coor, coor_tag])
        self.coor =  np.sort(self.coor, order=('y', 'x'))
        self.sites = len(self.coor['tag'])

    def plt_lattice(self, ms=30, fs=20, plt_label=False, figsize=None):
        '''
        Plot lattice.

        :param ms: Default value 30. Markersize. 
        :param lw: Default value 5. Linewidth of the hoppings.
        :param plt_label: Default value False. Plot site labels.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        colors = ['b', 'r', 'g', 'y', 'm', 'k']  # color plot 
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
        if plt_label:
            labels = ['{}'.format(i) for i in range(self.sites)]
            for l, x, y in zip(labels, self.coor['x'], self.coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        #lgnd = plt.legend(loc='upper right', numpoints=1)
        #for i, t in enumerate(self.tags): 
        #    lgnd.legendHandles[i]._legmarker.set_markersize(15)
        #plt.draw()
        return fig
