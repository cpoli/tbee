import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import os

class plotTB:
    '''
    Plot the results of class **eigTB**.

    :param sys: class instance **eigTB**.
    '''
    def __init__(self, sys):
        self.sys = sys
        self.colors = ['b', 'r', 'g', 'y', 'm', 'k']  # color plot 

    def plt_lattice(self, ms=30, lw=5, fs=20, c=3., plt_hop=False, plt_label=None, figsize=None):
        '''
        Plot lattice.

        :param hop: Default value []. Hoppings given by the class **eigTB**. 
            If not empty, plot the hoppings.
        :param ms: Default value 30. Markersize. 
        :param lw: Default value 5. Linewidth of the hoppings.
        :param fs: Default value 20. Fontsize.
        :param c: Default value 3. Coefficient, Bonds linewidth given by c*hop.
        :param label: Default value False. Plot site labels.
        :param figsize: Default value False. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots(figsize=figsize)
        # bonds
        for i in range(len(self.sys.hop)): 
            plt.plot([self.sys.coor['x'][self.sys.hop['i'][i]], self.sys.coor['x'][self.sys.hop['j'][i]]],
                        [self.sys.coor['y'][self.sys.hop['i'][i]], self.sys.coor['y'][self.sys.hop['j'][i]]],
                        'k', lw=c*self.sys.hop['t'][i].real)
        # sites
        for c, t in zip(self.colors, self.sys.tags):
            plt.plot(self.sys.coor['x'][self.sys.coor['tag'] == t],
                    self.sys.coor['y'][self.sys.coor['tag'] == t],
                    'o', color=c, ms=ms, markeredgecolor='none')
        ax.axis('off')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(self.sys.coor['x'])-0.5, np.max(self.sys.coor['x'])+0.5])
        ax.set_ylim([np.min(self.sys.coor['y'])-0.5, np.max(self.sys.coor['y'])+0.5])
        # labels
        if plt_label:
            labels = ['{}'.format(i) for i in range(self.sys.sites)]
            for l, x, y in zip(labels, self.sys.coor['x'], self.sys.coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        plt.draw()
        return fig

    def plt_spec(self, en_lims=[], pola=[], pola_tag='', ms=10, fs=20):
        '''
        Plot spectrum (eigenenergies real part (blue circles), 
        and sublattice polarization if *pola* not empty (red circles).

        :param en: Eigenenergies.
        :param pola: Defaut value []. Output of the class **eigTB**.
        :param pola_tag: Tag of the sublattice.
        :param ms: Default value 10. Markersize.
        :param fs: Default value 20. Fontsize.

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax1 = plt.subplots()
        x = np.arange(-self.sys.sites//2, self.sys.sites-self.sys.sites//2) + 1
        if en_lims == []:
            en_max = np.max(self.sys.en.real)
            ax1.set_ylim([-en_max-0.1, en_max+0.1])
            ind = np.ones(self.sys.sites, bool)
        else:
            ind = (self.sys.en > en_lims[0]) & (self.sys.en < en_lims[1])
            ax1.set_ylim([en_lims[0]-0.1, en_lims[1]+0.1]) 
        ax1.plot(x[ind], self.sys.en.real[ind], 'ob', markersize=ms)
        ax1.set_title('Spectrum', fontsize=fs)
        ax1.set_xlabel('$n$', fontsize=fs)
        ax1.set_ylabel('$E_n$', fontsize=fs, color='blue')
        for label in ax1.get_yticklabels():
            label.set_color('b')
        if pola_tag:
            i_tag = np.argwhere(self.sys.tags == pola_tag)
            ax2 = ax1.twinx()
            ax2.plot(x[ind], np.ravel(self.sys.pola[ind, i_tag]), 'or', markersize=(4*ms)//5)
            str_tag = pola_tag.decode('ascii')
            ylabel = '$<' + str_tag.upper() + '|' + str_tag.upper() + '>$' 
            ax2.set_ylabel(ylabel, fontsize=fs, color='red')
            ax2.set_ylim([-0.1, 1.1])
            ax2.set_yticks([0, 0.5, 1])
            if en_lims == []:
                ax2.set_xlim([x[0]-0.5, x[-1]+0.5])
            for label in ax2.get_yticklabels():
                label.set_color('r')
        plt.xlim([x[ind][0]-.1, x[ind][-1]+0.1])
        xa = ax1.get_xaxis()
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        plt.draw()
        return fig

    def plt_intensity1d(self, intensity, ms=20, lw=2, fs=20, title=''):
        '''
        Plot intensity.

        :param intensity: Intensity.
        :param ms: Default value 20. Markersize.
        :param lw: Default value 2. Linewith.
        :param fs: Default value 20. Font size.
        :param title: Default value ''. Figure title.
        '''
        fig, ax = plt.subplots()
        ax.set_xlabel('$j$', fontsize=fs)
        ax.set_ylabel('$|\psi^{(j)}|^2$', fontsize=fs)
        ax.set_title('Intensity', fontsize=fs)
        for t, c in zip(self.sys.tags, self.colors):
            plt.plot(self.sys.coor['x'][self.sys.coor['tag'] == t],
                        intensity[self.sys.coor['tag'] == t],
                        '-o', color=c, ms=ms, lw=lw)
        plt.xlim([-1., self.sys.sites])
        plt.ylim([0., np.max(intensity)+.02])
        plt.draw()
        return fig

    def plt_intensity(self, intensity, s=300, fs=20, title='$|\psi_n|^2$', lims=[]):
        '''
        Plot the intensity. Colormap with identical disk shape.

        :param intensity: Intensity.
        :param s: Default value 300. Disk size.
        :param fs: Default value 20. Font size.
        :param title: Default value '$|\psi_n|^2$'. Title.
        :param lims: Colormap limits. 

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots()
        plt.title(title, fontsize=fs)
        map_red = plt.get_cmap('Reds')
        if lims == []:
            lims = [0., np.max(intensity)]
            y_ticks = ['0', 'max']
        else: 
            y_ticks = lims
        plt.scatter(self.sys.coor['x'], self.sys.coor['y'], c=intensity, s=s, 
                         cmap=map_red, vmin=lims[0], vmax=lims[1])
        cbar = plt.colorbar(ticks=lims)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(np.min(self.sys.coor['x'])-1., np.max(self.sys.coor['x'])+1.)
        ax.set_ylim(np.min(self.sys.coor['y'])-1., np.max(self.sys.coor['y'])+1.)
        cbar.ax.set_yticklabels([y_ticks[0], y_ticks[1]])
        ax.set_aspect('equal')
        plt.draw()
        return fig

    def plt_intensity_disk(self, intensity, hop=[], s=1000., lw=1, fs=20, 
                                       title='$|\psi_n|^2$'):
        '''
        Plot the intensity. Intensity propotional to disk shape.

        :param intensity: Intensity.
        :param hop: Default value []. Hoppings. 
        :param s: Default value 1000. Circle size given by *s * intensity*.
        :param lw: Default value 1. Line width of the bonds.
        :param fs: Default value 20. Fontsize.
        :param title: Default value ''. Figure title.

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots()
        ax.set_xlabel('$i$', fontsize=fs)
        ax.set_ylabel('$j$', fontsize=fs)
        ax.set_title('$|\psi_n|^2$' + ' ' + add, fontsize=fs)
        if hop != []:
            for i in range(len(hop)): 
                plt.plot([self.coor['x'][hop['i'][i]], self.coor['x'][hop['j'][i]]],
                            [self.coor['y'][hop['i'][i]], self.coor['y'][hop['j'][i]]],
                            self.colors[0], lw=lw, c='k')

        for t, c in zip(self.tags, self.colors):
            plt.scatter(self.coor['x'][self.coor['tag'] == t],
                        self.coor['y'][self.coor['tag'] == t],
                        s=s*intensity[self.coor['tag'] == t],
                        c=c, alpha=0.5)

        ax.set_aspect('equal')
        ax.axis('off')
        ax.set_ylim([np.min(self.coor['y'])-1.5, np.max(self.coor['y'])+1.5])
        plt.draw()
        return fig


class saveFigTB():
    '''
    Create folder and save figures / animations obtained via **plotTB** or **propagationTB**.
    '''
    def __init__(self, sys, hop, hop_tags, dir_name, ext='png'):
        self.sys = sys
        self.hop_tags = hop_tags  # hopping tags
        self.hop = hop   # hoppings
        self.ext = ext  # file extension
        self.dir_main = '../TBfig/'  # main directory name
        self.dir_name = self.dir_name(dir_name)  # directory name
        self.check_dir()

    def dir_name(self, dir_name):
        '''
        Set the name of the directory in which the figures 
        are stored.

        :param dir_name: String. First part of the directory name. 
        '''
        dir_name = self.dir_main + dir_name + '_n{}'.format(self.sys.sites)
        return dir_name

    def check_dir(self):
        '''
        Create the directory to store the figures exists.
        '''
        if not os.path.exists(self.dir_main):
            os.makedirs(self.dir_main)
        if not os.path.exists(self.dir_name):
            os.makedirs(self.dir_name)

    def file_name(self):
        '''
        Create the file name.

        :returns:
            * **file_name** -- File name.
        '''
        file_name = ''
        for t, h in zip(self.hop_tags, self.hop): 
            file_name += '_' + t + str(complex(h+0)).replace('.', ',')
        for t, o in zip(self.sys.tags, self.sys.on): 
            file_name += '_e' + str(t)[2] + str(complex(o+0)).replace('.', ',')
        for key, val in self.sys.get_params().items():
            if val != 0.:
                file_name += '_' + key + '_' + str(val).replace('.', ',')
        return file_name

    def save_fig(self, fig, name=''):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  Fist part of the file name.
        '''
        name_file = self.dir_name + '/' + name + self.file_name() + '.' + self.ext
        fig.savefig(name_file, format=self.ext)
