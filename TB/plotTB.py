import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import os

def test_plot_lattice(tags, ms, fs, colors, plt_index, figsize):
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
            raise ValueError('\n\nParameter figsize must be list/tuple of length 2\n')
        if isinstance(figsize, (list, tuple)):
            if any([val < 0 for val in figsize]):
                raise ValueError('\n\nParameter figsize must contain\
                                        positive numbers.\n')


    return colors, figsize

def test_colors(colors, unit_cell):
    ''' 
    Test *colors* parameter.
    '''

    """
    if colors is None:
        pass
    else:
        if len(colors) != len(self.sys.lat.unit_cell):
            raise ValueError('\n\nlength of colors must be a least equal to\
                                        length of tags.\n')
        plt_colors = ['b', 'r', 'g', 'y', 'm', 'k', 'c']
      
        for c in colors:
            if not (c[0] == '#' and len(c) == 7) and c not in plt_colors:
                raise ValueError('\n\nParameter colors must be a list/tuple of\
                                          of hexagonal colors or matplotlib colors:\
                                          ["b", "g", "r", "c", "m", "y", "k"]\n')
    """                                   
class plotTB:
    '''
    Plot the results of class **eigTB**.

    :param sys: class instance **eigTB**.
    :param colors: Default value []. Color plot.
    '''
    def __init__(self, sys, colors=None):
        self.sys = sys
        if colors is None:
            self.colors = ['b', 'r', 'g', 'y', 'm', 'k']
        else:
            self.colors = colors
        test_colors(self.colors, self.sys.lat.unit_cell)

    def lattice_generic(self, fig, ax, coor, ms, lw, fs, c, plt_hop, plt_index, figsize):
        '''
        Private method.
        '''
        # plot sites
        for color, dic in zip(self.colors, self.sys.lat.unit_cell):
            plt.plot(coor['x'][coor['tag'] == dic['tag']],
                       coor['y'][coor['tag'] == dic['tag']],
                       'o', color=color, ms=ms, markeredgecolor='none')
        ax.axis('off')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(coor['x'])-0.5, np.max(coor['x'])+0.5])
        ax.set_ylim([np.min(coor['y'])-0.5, np.max(coor['y'])+0.5])
        # plot indices
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sys.lat.sites)]
            for l, x, y in zip(indices, coor['x'], coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                            textcoords='offset points', ha='right',
                            va='bottom', size=fs)
        plt.draw()
        return fig

    def lattice(self, ms=20, lw=5, fs=20, c=3., plt_hop=False, plt_index=None, figsize=None):
        '''
        Plot lattice.

        :param plt_hop: Default value None. Hoppings given by the class **eigTB**. 
            If not empty, plot the hoppings.
        :param ms: Default value 20. Markersize. 
        :param c: Default value 3. Coefficient. Bond linewidths given by c*hop.
        :param fs: Default value 20. Fontsize.
        :param plt_index: Default value False. Plot site labels.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots(figsize=figsize)
        # plot bonds
        if plt_hop:
            for i in range(len(self.sys.hop)):
                plt.plot([self.sys.lat.coor['x'][self.sys.hop['i'][i]], 
                             self.sys.lat.coor['x'][self.sys.hop['j'][i]]],
                            [self.sys.lat.coor['y'][self.sys.hop['i'][i]], 
                             self.sys.lat.coor['y'][self.sys.hop['j'][i]]],
                            'k', lw=c*self.sys.hop['t'][i].real) 
        return self.lattice_generic(fig, ax, self.sys.lat.coor, ms, lw, fs, c, plt_hop, plt_index, figsize)

    def lattice_hop(self, ms=20, lw=5, fs=20, c=3., plt_hop=False, plt_index=None, figsize=None):
        '''
        Plot lattice.

        :param plt_hop: Default value None. Hoppings given by the class **eigTB**. 
            If not empty, plot the hoppings.
        :param ms: Default value 20. Markersize. 
        :param c: Default value 3. Coefficient. Bond linewidths given by c*hop.
        :param fs: Default value 20. Fontsize.
        :param plt_index: Default value False. Plot site labels.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots(figsize=figsize)
        # bonds
        if plt_hop:
            for i in range(len(self.sys.hop)):
                plt.plot([self.sys.coor_hop['x'][self.sys.hop['i'][i]], self.sys.coor_hop['x'][self.sys.hop['j'][i]]],
                           [self.sys.coor_hop['y'][self.sys.hop['i'][i]], self.sys.coor_hop['y'][self.sys.hop['j'][i]]],
                           'k', lw=lw)
        return self.lattice_generic(fig, ax, self.sys.coor_hop, ms, lw, fs, c, plt_hop, plt_index, figsize) 

    def spectrum(self, en_lims=None, tag_pola=b'', ms=10, fs=20):
        '''
        Plot spectrum (eigenenergies real part (blue circles), 
        and sublattice polarization if *pola* not empty (red circles).

        :param en: np.array. Eigenenergies.
        :param en_lims: Default value []. Energy limits (size 2).
        :param tag_pola: Default value []. Byte type. Tag of the sublattice.
        :param ms: Default value 10. Markersize.
        :param fs: Default value 20. Fontsize.

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax1 = plt.subplots()
        ax1 = plt.gca()
        x_min = - np.floor(self.sys.lat.sites / 2.)
        x_max = np.ceil(self.sys.lat.sites / 2.)
        x = np.arange(x_min, x_max) 
        if en_lims == None:
            en_max = np.max(self.sys.en.real)
            ax1.set_ylim([-en_max-0.2, en_max+0.2])
            ind = np.ones(self.sys.lat.sites, bool)
        else:
            ind = (self.sys.en > en_lims[0]) & (self.sys.en < en_lims[1])
            ax1.set_ylim([en_lims[0]-0.1, en_lims[1]+0.1])
        ax1.plot(x[ind], self.sys.en.real[ind], 'ob', markersize=ms)
        ax1.set_title('Spectrum', fontsize=fs)
        ax1.set_xlabel('$n$', fontsize=fs)
        ax1.set_ylabel('$E_n$', fontsize=fs, color='blue')
        for label in ax1.get_yticklabels():
            label.set_color('b')
        if tag_pola:
            tags = [dic['tag'] for dic in self.sys.lat.unit_cell]
            i_tag= [ind for ind, val in enumerate(tags) if val == tag_pola]
            ax2 = plt.twinx()
            ax2.plot(x[ind], np.ravel(self.sys.pola[ind, i_tag]), 'or', markersize=(4*ms)//5)
            str_tag = tag_pola.decode('ascii')
            ylabel = '$<' + str_tag.upper() + '|' + str_tag.upper() + '>$' 
            ax2.set_ylabel(ylabel, fontsize=fs, color='red')
            ax2.set_ylim([-0.1, 1.1])
            ax2.set_yticks([0, 0.5, 1])
            for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_fontsize(18) 
            if en_lims == []:
                ax2.set_xlim([x[0]-0.5, x[-1]+0.5])
            for label in ax2.get_yticklabels():
                label.set_color('r')
        for label in ax1.xaxis.get_majorticklabels():
            label.set_fontsize(18) 
        for label in ax1.yaxis.get_majorticklabels():
            label.set_fontsize(18)
        xa = ax1.get_xaxis()
        ax1.set_xlim([x[ind][0]-0.1, x[ind][-1]+0.1]) 
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def spectrum_hist(self, nbr_bins=61, en_lims=[-1.5, 1.5]):
        """
        Plot the spectrum.
            
        :param nbr_bins: Default value 101. Number of bins of the histogram.
        :param ener_lim:  Default value [-1.5, 1.5]. list of the energy min and max.
        :param ms: Default value 10. Size of the markers.
        :param save: Default value False. Save the two figures.
        """
        ind_en = np.argwhere((self.sys.en > en_lims[0]) & (self.sys.en < en_lims[1]))
        ind_en = np.ravel(ind_en)
        en = self.sys.en[ind_en]
        fig, ax = plt.subplots()
        fig.canvas.set_window_title('Spectrum')
        n, bins, patches = plt.hist(en, bins=nbr_bins, color='#00008B')
        ax.set_xlabel('$E$', fontsize=20)
        ax.set_ylabel('number of states', fontsize=20)
        ax.set_ylim([0, np.max(n)])
        ax.set_xlim(en_lims)

    def intensity1d(self, intensity, ms=20, lw=2, fs=20, title='Intensity'):
        '''
        Plot intensity.

        :param intensity: np.array. Field intensity.
        :param ms: Default value 20. Markersize.
        :param lw: Default value 2. Linewith, connect sublattice sites.
        :param fs: Default value 20. Font size.
        :param title: Default value 'Intensity'. Figure title.
        '''
        fig, ax = plt.subplots()
        ax.set_xlabel('$j$', fontsize=fs)
        ax.set_ylabel('$|\psi^{(j)}|^2$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        for t, c in zip(self.sys.lat.tags, self.colors):
            plt.plot(self.sys.lat.coor['x'][self.sys.lat.coor['tag'] == t],
                        intensity[self.sys.lat.coor['tag'] == t],
                        '-o', color=c, ms=ms, lw=lw)
        plt.xlim([-1., self.sys.lat.sites])
        plt.ylim([0., np.max(intensity)+.02])
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def intensity_disk(self, intensity, s=200, fs=20, title=r'$|\psi|^2$', lims=[]):
        '''
        Plot the intensity. Colormap with identical disk shape.

        :param intensity: np.array.Field intensity.
        :param s: Default value 300. Disk size.
        :param fs: Default value 20. Font size.
        :param title: Default value '$|\psi_n|^2$'. Title.
        :param lims: Colormap limits. 

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots()
        plt.title(title, fontsize=fs+5)
        map_red = plt.get_cmap('Reds')
        if lims == []:
            lims = [0., np.max(intensity)]
            y_ticks = ['0', 'max']
        else: 
            y_ticks = lims
        plt.scatter(self.sys.lat.coor['x'], self.sys.lat.coor['y'], c=intensity, s=s, 
                         cmap=map_red, vmin=lims[0], vmax=lims[1])
        cbar = plt.colorbar(ticks=lims)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(np.min(self.sys.lat.coor['x'])-1., np.max(self.sys.lat.coor['x'])+1.)
        ax.set_ylim(np.min(self.sys.lat.coor['y'])-1., np.max(self.sys.lat.coor['y'])+1.)
        cbar.ax.set_yticklabels([y_ticks[0], y_ticks[1]])
        cbar.ax.tick_params(labelsize=fs)
        ax.set_aspect('equal')
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def intensity_area(self, intensity, plt_hop=None, s=1000., lw=1, fs=20, title=''):
        '''
        Plot the intensity. Intensity propotional to disk shape.

        :param intensity: np.array. Intensity.
        :param hop: Default value []. Hoppings. 
        :param s: Default value 1000. Circle size given by *s * intensity*.
        :param lw: Default value 1. Bond Line widths.
        :param fs: Default value 20. Fontsize.
        :param title: Default value '$|\psi_{ij}|^2$'. Figure title.

        :returns:
            * **fig** -- Figure.
        '''
        fig, ax = plt.subplots()
        ax.set_xlabel('$i$', fontsize=fs)
        ax.set_ylabel('$j$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        if plt_hop:
            for i in range(len(self.sys.hop)):
                plt.plot([self.sys.lat.coor['x'][self.sys.hop['i'][i]], self.sys.lat.coor['x'][self.sys.hop['j'][i]]],
                           [self.sys.lat.coor['y'][self.sys.hop['i'][i]], self.sys.lat.coor['y'][self.sys.hop['j'][i]]],
                            'k', lw=lw)
        tags = [dic['tag'] for dic in self.sys.lat.unit_cell]
        for tag, color in zip(tags, self.colors):
            plt.scatter(self.sys.lat.coor['x'][self.sys.lat.coor['tag'] == tag],
                        self.sys.lat.coor['y'][self.sys.lat.coor['tag'] == tag],
                        s=100*s*intensity[self.sys.lat.coor['tag'] == tag],
                        c=color, alpha=0.5)
        ax.set_aspect('equal')
        ax.axis('off')
        x_lim = [np.min(self.sys.lat.coor['x'])-2., np.max(self.sys.lat.coor['x'])+2.]
        y_lim = [np.min(self.sys.lat.coor['y'])-2., np.max(self.sys.lat.coor['y'])+2.]
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def butterfly(self, en_lims=[-2., 2.], fs=20):
        '''
        Plot energies depending on strain.
        '''
        i_beta_min = np.argmin(np.abs(self.sys.betas))
        ind_en = np.argwhere((self.sys.butterfly[i_beta_min, :] > en_lims[0]) & 
                                        (self.sys.butterfly[i_beta_min, :] < en_lims[1]))
        ind_en = np.ravel(ind_en)
        fig, ax = plt.subplots()
        plt.title('Energies depending on strain', fontsize=fs)
        plt.xlabel(r'$\beta/\beta_{max}$', fontsize=fs)
        plt.ylabel('$E$', fontsize=fs)
        plt.yticks(np.arange(en_lims[0], en_lims[1]+1e-2), fontsize=fs)
        beta_max = max(self.sys.betas)
        plt.xticks((-beta_max, -0.5*beta_max, 0, 0.5*beta_max, beta_max), 
                        fontsize=fs)
        ax.set_xticklabels(('-1', '-1/2', '0', '1/2', '1'))
        plt.xlim([self.sys.betas[0], self.sys.betas[-1]])
        plt.ylim(en_lims)
        no_en = len(ind_en)
        for i in ind_en:
            plt.plot(self.sys.betas, self.sys.butterfly[:, i], 'b')
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def show(self):
        """
        Emulate Matplotlib method plt.show().
        """
        plt.show()


class saveTB():
    '''
    Create folder and save figures / animations obtained via 
    **latticeTB** , **plotTB** or **propagationTB**.
    '''
    def __init__(self, sys, dir_name, params={}, ext='png'):
        self.sys = sys
        self.params = params  # parameters dictionary
        self.ext = ext  # file extension
        self.dir_main = '../TBfig/'  # main directory name
        self.dir_name = self.dir_name(dir_name)  # directory name
        self.check_dir()

    def dir_name(self, dir_name):
        '''
        Set the name of the directory in which the figures are stored.

        :param dir_name: String. First part of the directory name. 
        '''
        dir_name = self.dir_main + dir_name + '_n{}'.format(self.sys.lat.sites)
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
        tags = [dic['tag'] for dic in self.sys.lat.unit_cell]
        for t, o in zip(tags, self.sys.ons): 
            file_name += '_e' + str(t)[2] + str(complex(o+0)).replace('.', ',')
        for key, val in self.params.items():
            file_name += '_' + key + str(complex(val+0)).replace('.', ',')
        return file_name

    def fig(self, fig, name=''):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. Fist part of the file name.
        '''
        name_file = self.dir_name + '/' + name + self.file_name() + '.' + self.ext
        fig.savefig(name_file, format=self.ext)

    def fig_lat(self, fig, name=''):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. First part of the file name.
        '''
        name_file = self.dir_name + '/' + name + '.' + self.ext
        fig.savefig(name_file, format=self.ext)

    def ani(self, ani, name='', fps=10):
        name_file = self.dir_name + '/' + name + '.mp4'
        ani.save(name_file, fps=fps, extra_args=['-vcodec', 'libx264'])
