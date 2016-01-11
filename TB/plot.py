import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import error_handling
import os


#################################
# CLASS PLOT
#################################


class plot:
    def __init__(self, sys, colors=None):
        '''
        Plot the results of the classes **lattice** or **system**.

        :param sys: class instance **system**.
        :param colors: Default value None. Color plot.
        '''
        error_handling.sys(sys) 
        self.sys = sys
        if colors is None:
            self.colors = ['b', 'r', 'g', 'y', 'm', 'k']
        else:
            self.colors = colors

    def plt_hopping(self, coor, hop, c):
        for i in range(len(hop)):
            plt.plot([coor['x'][hop['i'][i]], 
                        coor['x'][hop['j'][i]]],
                        [coor['y'][hop['i'][i]], 
                         coor['y'][hop['j'][i]]],
                        'k', lw=c*hop['t'][i].real)        

    def lattice_generic(self, coor, ms, lw, fs, c, plt_hop, plt_hop_low, 
                                    plt_index, axis, figsize):
        '''
        Private method.
        '''
        error_handling.positive_int(ms, 'ms')
        error_handling.positive_real(lw, 'lw')
        error_handling.positive_int(fs, 'fs')
        error_handling.positive_real(c, 'c')
        error_handling.boolean(plt_hop, 'plt_hop')
        error_handling.boolean(plt_hop_low, 'plt_hop_low')
        error_handling.boolean(plt_index, 'plt_index')
        error_handling.boolean(axis, 'axis')
        if figsize is None:
            figsize = (5, 4)
        error_handling.list_tuple_2elem(figsize, 'figsize')
        error_handling.positive_int(figsize[0], 'figsize[0]')
        error_handling.positive_int(figsize[1], 'figsize[1]')
        fig, ax = plt.subplots(figsize=figsize)
        # hoppings
        if plt_hop:
            self.plt_hopping(coor, self.sys.hop, c)
        if plt_hop_low:
            self.plt_hopping(coor, self.sys.hop_low, c)

        # plot sites
        for color, tag in zip(self.colors, self.sys.lat.tags):
            plt.plot(coor['x'][coor['tag'] == tag],
                       coor['y'][coor['tag'] == tag],
                       'o', color=color, ms=ms, markeredgecolor='none')
        ax.set_aspect('equal')
        ax.set_xlim([np.min(coor['x'])-1., np.max(coor['x'])+1.])
        ax.set_ylim([np.min(coor['y'])-1., np.max(coor['y'])+1.])
        if not axis:
            ax.axis('off')
        # plot indices
        if plt_index:
            indices = ['{}'.format(i) for i in range(self.sys.lat.sites)]
            for l, x, y in zip(indices, coor['x'], coor['y']):
                plt.annotate(l, xy=(x, y), xytext=(0, 0),
                                    textcoords='offset points',
                                    ha='right', va='bottom', size=fs)
        plt.draw()
        return fig

    def lattice(self, ms=20, lw=5, fs=20, c=3., plt_hop=False, plt_hop_low=False,
                    plt_index=False, axis=False, figsize=None):
        '''
        Plot lattice.

        :param ms: Default value 20. Markersize. 
        :param c: Default value 3. Coefficient. Hopping linewidths given by c*hop['t'].
        :param fs: Default value 20. Fontsize.
        :param plt_hop: Default value False. Plot hoppings.
        :param plt_hop_low: Default value False. Plot hoppings diagonal low.
        :param plt_index: Default value False. Plot site labels.
        :param axis: Default value False. Plot axis.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_coor(self.sys.lat.coor)
        return self.lattice_generic(self.sys.lat.coor, ms, lw, fs, c, plt_hop,
                                                 plt_hop_low, plt_index, axis, figsize)


    def lattice_hop(self, ms=20, lw=5, fs=20, c=3., plt_hop=False, plt_hop_low=False,
                            plt_index=False, axis=False, figsize=None):
        '''
        Plot lattice in hopping space.

        :param ms: Default value 20. Markersize. 
        :param c: Default value 3. Coefficient. Hopping linewidths given by c*hop['t'].
        :param fs: Default value 20. Fontsize.
        :param plt_hop: Default value False. Plot hoppings.
        :param plt_index: Default value False. Plot site labels.
        :param axis: Default value False. Plot axis.
        :param figsize: Default value None. Figsize. 

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_coor_hop(self.sys.coor_hop)
        return self.lattice_generic(self.sys.coor_hop, ms, lw, fs, c, plt_hop, 
                                                 plt_hop_low, plt_index, axis, figsize) 

    def spectrum(self, ms=10, fs=20, lims=None, tag_pola=None, ipr=None):
        '''
        Plot spectrum (eigenenergies real part (blue circles), 
        and sublattice polarization if *pola* not empty (red circles).

        :param ms: Default value 10. Markersize.
        :param fs: Default value 20. Fontsize.
        :param lims: List, lims[0] energy min, lims[1] energy max.
        :param tag_pola: Default value None. Binary char. Tag of the sublattice.
        :param ipr: Default value None. np.ndarray. Inverse Partitipation Ration.
  

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_en(self.sys.en)
        error_handling.positive_int(ms, 'ms')
        error_handling.positive_int(fs, 'fs')
        error_handling.lims(lims)
        fig, ax1 = plt.subplots()
        ax1 = plt.gca()
        x = np.arange(self.sys.lat.sites) 
        if lims is None:
            en_max = np.max(self.sys.en.real)
            ax1.set_ylim([-en_max-0.2, en_max+0.2])
            ind = np.ones(self.sys.lat.sites, bool)
        else:
            ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
            ax1.set_ylim([lims[0]-0.1, lims[1]+0.1])
        ax1.plot(x[ind], self.sys.en.real[ind], 'ob', markersize=ms)
        ax1.set_title('Spectrum', fontsize=fs)
        ax1.set_xlabel('$n$', fontsize=fs)
        ax1.set_ylabel('$E_n$', fontsize=fs, color='blue')
        for label in ax1.get_yticklabels():
            label.set_color('b')
        if tag_pola:
            error_handling.tag(tag_pola, self.sys.lat.tags)
            fig, ax2 = self.plt_pola(fig=fig, ax1=ax1, ms=ms, fs=fs, tag_pola=tag_pola, ind=ind)
        elif ipr:
            fig, ax2 = self.plt_ipr(fig=fig, ax1=ax1, ms=ms, fs=fs, tag_pola=tag_pola, ind=ind)
        for label in ax1.xaxis.get_majorticklabels():
            label.set_fontsize(fs) 
        for label in ax1.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        xa = ax1.get_xaxis()
        ax1.set_xlim([x[ind][0]-0.1, x[ind][-1]+0.1]) 
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def plt_pola(self, fig=None, ax1=None, ms=10, fs=20, lims=None, tag_pola=None, ind=None):
        if fig is None:
            error_handling.sys(self.sys)
            error_handling.empty_en(self.sys.en)
            error_handling.positive_int(ms, 'ms')
            error_handling.positive_int(fs, 'fs')
            error_handling.lims(lims)
            fig, ax2 = plt.subplots()
            ax2 = plt.gca()
            if lims is None:
                en_max = np.max(self.sys.en.real)
                ax2.set_ylim([-0.1, 1.1])
                ind = np.ones(self.sys.lat.sites, bool)
            else:
                ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
                ax2.set_ylim([lims[0]-0.1, lims[1]+0.1])
        else:
            ax2 = plt.twinx()
        error_handling.empty_pola(self.sys.pola)
        error_handling.tag(tag_pola, self.sys.lat.tags)
        x = np.arange(self.sys.lat.sites) 
        i_tag = self.sys.lat.tags == tag_pola
        ax2.plot(x[ind], np.ravel(self.sys.pola[ind, i_tag]), 'or', markersize=(4*ms)//5)
        str_tag = tag_pola.decode('ascii')
        ylabel = '$<' + str_tag.upper() + '|' + str_tag.upper() + '>$' 
        ax2.set_ylabel(ylabel, fontsize=fs, color='red')
        #ax2.set_ylim([-0.1, 1.1])
        ax2.set_yticks([0, 0.5, 1])
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(fs) 
        for label in ax2.get_yticklabels():
            label.set_color('r')
        return fig, ax2

    def plt_ipr(self, fig=None, ax1=None, ms=10, fs=20, lims=None, tag_pola=None, ind=None):
        if fig is None:
            error_handling.sys(self.sys)
            error_handling.empty_en(self.sys.en)
            error_handling.positive_int(ms, 'ms')
            error_handling.positive_int(fs, 'fs')
            error_handling.lims(lims)
            fig, ax2 = plt.subplots()
            ax2 = plt.gca()
            if lims is None:
                en_max = np.max(self.sys.en.real)
                ax2.set_ylim([-0.1, 1.1])
                ind = np.ones(self.sys.lat.sites, bool)
            else:
                ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
                ax2.set_ylim([lims[0]-0.1, lims[1]+0.1])
        else:
            ax2 = plt.twinx()
        error_handling.empty_ipr(self.sys.ipr)
        x = np.arange(self.sys.lat.sites) 
        ax2.plot(x[ind], self.sys.ipr[ind], 'or', markersize=(4*ms)//5)
        ax2.set_ylabel( 'IPR' , fontsize=fs, color='red')
        for tick in ax2.xaxis.get_major_ticks():
            tick.label.set_fontsize(fs) 
        for label in ax2.get_yticklabels():
            label.set_color('r')
        return fig, ax2

    def spectrum_hist(self, nbr_bins=61, lims=None):
        """
        Plot the spectrum.
            
        :param nbr_bins: Default value 101. Number of bins of the histogram.
        :param ener_lim:  Default value [-1.5, 1.5]. list of the energy min and max.
        :param ms: Default value 10. Size of the markers.
        :param save: Default value False. Save the two figures.
        """
        error_handling.empty_en(self.sys.en)
        error_handling.positive_int(nbr_bins, 'nbr_bins')
        error_handling.lims(lims)
        ind_en = np.argwhere((self.sys.en > lims[0]) & (self.sys.en < lims[1]))
        ind_en = np.ravel(ind_en)
        en = self.sys.en[ind_en]
        fig, ax = plt.subplots()
        fig.canvas.set_window_title('Spectrum')
        n, bins, patches = plt.hist(en, bins=nbr_bins, color='#00008B')
        ax.set_xlabel('$E$', fontsize=20)
        ax.set_ylabel('number of states', fontsize=20)
        ax.set_ylim([0, np.max(n)])
        ax.set_xlim(lims)

    def spectrum_complex(self, ms=10, fs=20, lims=None):
        '''
        Plot spectrum (eigenenergies real part (blue circles), 
        eigenenergies imaginary part (red circles).

        :param ms: Default value 10. Markersize.
        :param fs: Default value 20. Fontsize.
        :param lims: List, lims[0] energy min, lims[1] energy max.  

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_en(self.sys.en)
        error_handling.positive_int(ms, 'ms')
        error_handling.positive_int(fs, 'fs')
        error_handling.lims(lims)
        fig, ax1 = plt.subplots()
        ax1 = plt.gca()
        x = np.arange(self.sys.lat.sites) 
        if lims is None:
            en_max = np.max(self.sys.en.real)
            ax1.set_ylim([-en_max-0.2, en_max+0.2])
            ind = np.ones(self.sys.lat.sites, bool)
        else:
            ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
            ax1.set_ylim([lims[0]-0.1, lims[1]+0.1])
        ax1.plot(x[ind], self.sys.en.real[ind], 'ob', markersize=ms)
        ax1.plot(x[ind], self.sys.en.imag[ind], 'or', markersize=ms)
        ax1.set_title('Spectrum', fontsize=fs)
        ax1.set_xlabel('$n$', fontsize=fs)
        ax1.set_ylabel('$E_n$', fontsize=fs, color='blue')
        for label in ax1.xaxis.get_majorticklabels():
            label.set_fontsize(fs) 
        for label in ax1.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        xa = ax1.get_xaxis()
        ax1.set_xlim([x[ind][0]-0.1, x[ind][-1]+0.1]) 
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def intensity_1d(self, intensity, ms=20, lw=2, fs=20, title=r'$|\psi^{(j)}|^2$'):
        '''
        Plot intensity.

        :param intensity: np.array. Field intensity.
        :param ms: Default value 20. Markersize.
        :param lw: Default value 2. Linewith, connect sublattice sites.
        :param fs: Default value 20. Font size.
        :param title: Default value 'Intensity'. Figure title.
        '''
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.ndarray_null(intensity, 'intensity')
        error_handling.empty_coor(self.sys.lat.coor)
        error_handling.positive_int(ms, 'ms')
        error_handling.positive_int(lw, 'lw')
        error_handling.positive_int(fs, 'fs')
        error_handling.string(title, 'title')
        fig, ax = plt.subplots()
        ax.set_xlabel('$j$', fontsize=fs)
        ax.set_ylabel(title, fontsize=fs)
        ax.set_title(title, fontsize=fs)
        for t, c in zip(self.sys.lat.tags, self.colors):
            plt.plot(self.sys.lat.coor['x'][self.sys.lat.coor['tag'] == t],
                        intensity[self.sys.lat.coor['tag'] == t],
                        '-o', color=c, ms=ms, lw=lw)
        plt.xlim([-1., self.sys.lat.sites])
        plt.ylim([0., np.max(intensity)+.05])
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def intensity_disk(self, intensity, s=200, fs=20, lims=None, title=r'$|\psi|^2$'):
        '''
        Plot the intensity. Colormap with identical disk shape.

        :param intensity: np.array.Field intensity.
        :param s: Default value 300. Disk size.
        :param fs: Default value 20. Font size.
        :param lims: Colormap limits. 
        :param title: Default value '$|\psi_n|^2$'. Title.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.sys(self.sys)
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.ndarray_null(intensity, 'intensity')
        error_handling.empty_coor(self.sys.lat.coor)
        error_handling.positive_int(s, 's')
        error_handling.positive_int(fs, 'fs')
        error_handling.string(title, 'title')
        fig, ax = plt.subplots()
        plt.title(title, fontsize=fs+5)
        map_red = plt.get_cmap('Reds')
        if lims is None:
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

    def intensity_area(self, intensity, s=1000., lw=1, fs=20, plt_hop=False, title=r'$|\psi|^2$'):
        '''
        Plot the intensity. Intensity propotional to disk shape.

        :param intensity: np.array. Intensity.
        :param plt_hop: Default value False. Plot hoppings. 
        :param s: Default value 1000. Circle size given by :math:`s * intensity`.
        :param lw: Default value 1. Hopping linewidths.
        :param fs: Default value 20. Fontsize.
        :param title: Default value '$|\psi_{ij}|^2$'. Figure title.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.ndarray_null(intensity, 'intensity')
        error_handling.empty_coor(self.sys.lat.coor)
        error_handling.positive_int(s, 's')
        error_handling.positive_int(fs, 'fs')
        error_handling.boolean(plt_hop, 'plt_hop')
        error_handling.string(title, 'title')
        fig, ax = plt.subplots()
        ax.set_xlabel('$i$', fontsize=fs)
        ax.set_ylabel('$j$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        if plt_hop:
            for i in range(len(self.sys.hop)):
                plt.plot([self.sys.lat.coor['x'][self.sys.hop['i'][i]], 
                            self.sys.lat.coor['x'][self.sys.hop['j'][i]]],
                           [self.sys.lat.coor['y'][self.sys.hop['i'][i]], 
                            self.sys.lat.coor['y'][self.sys.hop['j'][i]]],
                            'k', lw=lw)
        for tag, color in zip(self.sys.lat.tags, self.colors):
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

    def butterfly(self, lw=1, fs=20, lims=[-2., 2.], title=''):
        '''
        Plot energies depending on a parameter.
        '''
        #if 'x' in locals() or 'x' in globals():
        #    raise NameError('butterfly is not defined')
        error_handling.positive_int(lw, 'lw')
        error_handling.positive_int(fs, 'fs')
        error_handling.string(title, 'title')

        i_beta_min = np.argmin(np.abs(self.sys.betas))
        ind_en = np.argwhere((self.sys.butterfly[i_beta_min, :] > lims[0]) & 
                                        (self.sys.butterfly[i_beta_min, :] < lims[1]))
        ind_en = np.ravel(ind_en)
        fig, ax = plt.subplots()
        plt.title('Energies depending on strain', fontsize=fs)
        plt.xlabel(r'$\beta/\beta_{max}$', fontsize=fs)
        plt.ylabel('$E$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        plt.yticks(np.arange(lims[0], lims[1]+1e-2), fontsize=fs)
        beta_max = max(self.sys.betas)
        plt.xticks([-beta_max, -0.5*beta_max, 0, 
                        0.5*beta_max, beta_max], fontsize=fs)
        ax.set_xticklabels(('-1', '-1/2', '0', '1/2', '1'))
        plt.xlim([self.sys.betas[0], self.sys.betas[-1]])
        plt.ylim(lims)
        no_en = len(ind_en)
        for i in ind_en:
            plt.plot(self.sys.betas, self.sys.butterfly[:, i], 'b', lw=lw)
        fig.set_tight_layout(True)
        plt.draw()
        return fig

    def show(self):
        """
        Emulate Matplotlib method plt.show().
        """
        plt.show()


#################################
# CLASS SAVE
#################################


class save():
    def __init__(self, dir_name, dir_main=None, params={}, file_format='png'):
        '''
        Create folder and save figures / animations obtained via 
        **plot** or **propagation**. 
        Plot the results of the classes **lattice** or **system**.

        :param dir_main: Name of the directory.
        :param dir_name: Default value None. Relative path of the main directory.
          if None, figures stored in ''../TBfig/'dir_name/'
        :param params: dictionary. file name information
        :param file_format: Default value 'png'. Figure format.
        '''
        error_handling.string(dir_name, 'dir_name')
        error_handling.string(dir_main, 'dir_main')
        error_handling.file_format(file_format)
        self.params = params  # dictionary parameters
        self.file_format = file_format  # file format
        if dir_main is None:
            self.dir_main = '../TBfig/'  # main directory name
        else:
            self.dir_main = dir_name
        self.dir_name = self.dir_name(dir_name)  # directory name
        self.check_dir()

    def dir_name(self, dir_name):
        '''
        Set the name of the directory in which the figures are stored.

        :param dir_name: String. First part of the directory name. 
        '''

        dir_name = self.dir_main + dir_name
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
        for key, val in self.params.items():
            file_name += '_' + key + str(complex(val+0)).replace('.', ',')
        return file_name

    def fig(self, fig, name):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. Fist part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + self.file_name() + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def fig_lat(self, fig, name):
        '''
        Save the figure in the directory defined by the method *dir_name()*.

        :param fig: Matplotlib fig.
        :param name:  String. First part of the file name.
        '''
        error_handling.fig(fig)
        error_handling.string(name, 'name')
        name_file = self.dir_name + '/' + name + '.' + self.file_format
        fig.savefig(name_file, format=self.file_format)

    def ani(self, ani, name, fps=10):
        error_handling.ani(ani)
        error_handling.string(name, 'name')
        error_handling.positive_int(fps, 'fps')
        name_file = self.dir_name + '/' + name + '.mp4'
        ani.save(name_file, fps=fps, extra_args=['-vcodec', 'libx264'])



"""
    def __init__(self, dir_name, dir_main=None, lat=None, sys=None, params={}, ext='png'):
        test_string(dir_name)
        test_string(dir_main)
        test_string(ext)
        self.params = params  # dictionary parameters
        self.lat = lat
        self.sys = sys
        if lat:
            self.sites = self.lat.sites
            self.tags = self.lat.tags
        else:
            self.sites = self.sys.lat.sites
            self.tags = self.sys.lat.tags
        self.ext = ext  # file extension
        if dir_main is None:
            self.dir_main = '../TBfig/'  # main directory name
        else:
            self.dir_main = dir_name
        self.dir_name = self.dir_name(dir_name)  # directory name
        self.check_dir()

    def dir_name(self, dir_name):
        '''
        Set the name of the directory in which the figures are stored.

        :param dir_name: String. First part of the directory name. 
        '''

        dir_name = self.dir_main + dir_name + '_n{}'.format(self.sites)
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
        if sys:
            for t, o in zip(self.tags, self.sys.ons): 
                file_name += '_e' + str(t)[2] + str(complex(o+0)).replace('.', ',')
            for key, val in self.params.items():
                file_name += '_' + key + str(complex(val+0)).replace('.', ',')
        return file_name


"""