from __future__ import annotations

import numpy as np
from numpy.typing import NDArray, ArrayLike
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
import tbee.error_handling as error_handling
import tbee.dos as dos
from tbee.system import System
import os


#################################
# CLASS PLOT
#################################


class Plot:
    '''
    Plot the results of the classes **lattice** or **system**.

    :param sys: class instance **system**.
    :param colors: Default value None. Color plot.
    '''

    def __init__(self, sys: System, colors: list[str] | None = None) -> None:
        error_handling.sys(sys)
        self.sys = sys
        if colors is None:
            self.colors = ['b', 'r', 'g', 'y', 'm', 'k']
        else:
            self.colors = colors

    def plt_hopping(self, coor: NDArray, hop: NDArray, c: float) -> None:
        '''
        Private method called by *lattice_generic*.
        '''
        for i in range(len(hop)):
            plt.plot([coor['x'][hop['i'][i]],
                        coor['x'][hop['j'][i]]],
                        [coor['y'][hop['i'][i]],
                         coor['y'][hop['j'][i]]],
                        'k', lw=c*hop['t'][i].real)

    def lattice_generic(
        self,
        coor: NDArray,
        ms: float,
        lw: float,
        c: float,
        fs: float,
        axis: bool,
        plt_hop: bool,
        plt_hop_low: bool,
        plt_index: bool,
        figsize: tuple[float, float] | None,
    ) -> Figure:
        '''
        Private method called by *lattice* and *lattice_hop*.
        '''
        error_handling.positive_real(ms, 'ms')
        error_handling.positive_real(lw, 'lw')
        error_handling.positive_real(c, 'c')
        error_handling.positive_real(fs, 'fs')
        error_handling.boolean(axis, 'axis')
        error_handling.boolean(plt_hop, 'plt_hop')
        error_handling.boolean(plt_hop_low, 'plt_hop_low')
        error_handling.boolean(plt_index, 'plt_index')
        error_handling.tuple_2elem(figsize, 'figsize')
        fig, ax = plt.subplots(figsize=figsize)
        # hoppings
        if plt_hop:
            error_handling.empty_ndarray(self.sys.hop, 'sys.hop')
            self.plt_hopping(coor, self.sys.hop[self.sys.hop['ang']>=0], c)
        if plt_hop_low:
            error_handling.empty_ndarray(self.sys.hop, 'sys.hop')
            self.plt_hopping(coor, self.sys.hop[self.sys.hop['ang']<0], c)
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

    def lattice(
        self,
        ms: float = 20,
        lw: float = 5.,
        c: float = 3.,
        fs: float = 20,
        axis: bool = False,
        plt_hop: bool = False,
        plt_hop_low: bool = False,
        plt_index: bool = False,
        figsize: tuple[float, float] | None = None,
    ) -> Figure:
        '''
        Plot lattice.

        :param ms: Positive number. Default value 20. Markersize.
        :param c: Positive number. Default value 3. 
            Coefficient. Hopping linewidths given by c*hop['t'].
        :param fs: Positive number. Default value 20. Fontsize.
        :param plt_hop: Boolean. Default value False. Plot hoppings.
        :param plt_hop_low: Boolean. Default value False. 
            Plot hoppings diagonal low.
        :param plt_index: Boolean. Default value False. Plot site labels.
        :param axis: Boolean. Default value False. Plot axis.
        :param figsize: Tuple. Default value None. Figure size.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.lat.coor, 'lat.get_lattice')
        return self.lattice_generic(self.sys.lat.coor, ms, lw, c, fs, axis, plt_hop,
                                                 plt_hop_low, plt_index, figsize)

    def lattice_hop(
        self,
        ms: float = 20,
        lw: float = 5,
        c: float = 3.,
        fs: float = 20,
        axis: bool = False,
        plt_hop: bool = False,
        plt_hop_low: bool = False,
        plt_index: bool = False,
        figsize: tuple[float, float] | None = None,
    ) -> Figure:
        '''
        Plot lattice in hopping space.

        :param ms: Positive Float. Default value 20. Markersize.
        :param c: Positive Float. Default value 3. Coefficient.
            Hopping linewidths given by c*hop['t'].
        :param fs: Positive Float. Default value 20. Fontsize.
        :param axis: Boolean. Default value False. Plot axis.
        :param plt_hop: Boolean. Default value False. Plot hoppings.
        :param plt_index: Boolean. Default value False. Plot site labels.
        :param figsize: Tuple. Default value None. Figure size.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.coor_hop, 'sys.get_coor_hop')
        return self.lattice_generic(self.sys.coor_hop, ms, lw, c, fs, axis, plt_hop,
                                                 plt_hop_low, plt_index, figsize)


    def spectrum_hist(self, nbr_bins: int = 61, fs: float = 20, lims: tuple[float, float] | None = None) -> None:
        """
        Plot the spectrum.
            
        :param nbr_bins: Default value 101. Number of bins of the histogram.
        :param lims: List, lims[0] energy min, lims[1] energy max.
        """
        error_handling.empty_ndarray(self.sys.en, 'sys.get_eig')
        error_handling.positive_real(nbr_bins, 'nbr_bins')
        error_handling.lims(lims)
        fig, ax = plt.subplots()
        if lims is None:
            en_max = np.max(self.sys.en.real)
            ind_en = np.ones(self.sys.lat.sites, bool)
            ax.set_ylim([-en_max, en_max])
        else:
            ind_en = np.argwhere((self.sys.en > lims[0]) & (self.sys.en < lims[1]))
            ind_en = np.ravel(ind_en)
            ax.set_xlim(lims)
        en = self.sys.en[ind_en]
        n, bins, patches = plt.hist(en, bins=nbr_bins, color='b', alpha=0.8)
        ax.set_title('Spectrum', fontsize=fs)
        ax.set_xlabel('$E$', fontsize=fs)
        ax.set_ylabel('number of states', fontsize=fs)
        ax.set_ylim([0, np.max(n)+1])
        for label in ax.xaxis.get_majorticklabels():
            label.set_fontsize(fs)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(fs)

    def dos(
        self,
        broadening: float = 0.05,
        kernel: str = 'gaussian',
        e_grid: ArrayLike | None = None,
        fs: float = 20,
        lw: float = 2.,
        figsize: tuple[float, float] | None = None,
    ) -> Figure:
        '''
        Plot the (broadened) density of states, see *tbee.dos.density_of_states*.

        :param broadening: Positive real number. Default value 0.05. Kernel width.
        :param kernel: String. Default value 'gaussian'. 'gaussian' or 'lorentzian'.
        :param e_grid: Real ndarray. Default value None. Energies at which to
            evaluate the density of states.
        :param fs: Positive number. Default value 20. Fontsize.
        :param lw: Positive number. Default value 2. Linewidth.
        :param figsize: Tuple. Default value None. Figure size.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.en, 'sys.get_eig')
        error_handling.positive_real(fs, 'fs')
        error_handling.positive_real(lw, 'lw')
        error_handling.tuple_2elem(figsize, 'figsize')
        e_grid, rho = dos.density_of_states(self.sys.en, e_grid=e_grid,
                                                                   broadening=broadening, kernel=kernel)
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(e_grid, rho, 'b', lw=lw)
        ax.fill_between(e_grid, rho, color='b', alpha=0.2)
        ax.set_xlim([e_grid[0], e_grid[-1]])
        ax.set_ylim([0., None])
        ax.set_title('Density of states', fontsize=fs)
        ax.set_xlabel('$E$', fontsize=fs)
        ax.set_ylabel(r'$\rho(E)$', fontsize=fs)
        for label in ax.xaxis.get_majorticklabels():
            label.set_fontsize(fs)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def spectrum(
        self,
        ms: float = 10,
        fs: float = 20,
        lims: tuple[float, float] | None = None,
        tag_pola: str | None = None,
        ipr: bool | None = None,
        peterman: bool | None = None,
    ) -> Figure:
        '''
        Plot spectrum (eigenenergies real part (blue circles),
        and sublattice polarization if *pola* not empty (red circles).

        :param ms: Default value 10. Markersize.
        :param fs: Default value 20. Fontsize.
        :param lims: List, lims[0] energy min, lims[1] energy max.
        :param tag_pola: Default value None. One-character string. Tag of the sublattice.
        :param ipr: Default value None. If True plot the Inverse Partitipation Ration.
        :param petermann: Default value None. If True plot the Petermann factor.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.en, 'sys.get_eig')
        error_handling.positive_real(ms, 'ms')
        error_handling.positive_real(fs, 'fs')
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
            fig, ax2 = self.polarization(fig=fig, ax1=ax1, ms=ms, fs=fs, tag_pola=tag_pola, ind=ind)
        elif ipr:
            fig, ax2 = self.ipr(fig=fig, ax1=ax1, ms=ms, fs=fs, ind=ind)
        elif peterman:
            fig, ax2 = self.petermann(fig=fig, ax1=ax1, ms=ms, fs=fs, ind=ind)
        for label in ax1.xaxis.get_majorticklabels():
            label.set_fontsize(fs)
        for label in ax1.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        xa = ax1.get_xaxis()
        ax1.set_xlim([x[ind][0]-0.5, x[ind][-1]+0.5])
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def polarization(
        self,
        fig: Figure | None = None,
        ax1: Axes | None = None,
        ms: float = 10.,
        fs: float = 20.,
        lims: tuple[float, float] | None = None,
        tag_pola: str | None = None,
        ind: NDArray | None = None,
    ) -> tuple[Figure, Axes]:
        '''
        Plot sublattice polarization.

        :param fig: Figure. Default value None. (used by the method spectrum).
        :param ax1: Axis. Default value None. (used by the method spectrum).
        :param ms: Positive Float. Default value 10. Markersize.
        :param fs: Positive Float. Default value 20. Fontsize.
        :param lims: List, lims[0] energy min, lims[1] energy max.
        :param tag_pola: One-character string. Default value None. Tag of the sublattice.
        :param ind: List. Default value None. List of indices. (used in the method spectrum).

        :returns:
            * **fig** -- Figure.
        '''
        if fig is None:
            error_handling.sys(self.sys)
            error_handling.empty_ndarray(self.sys.en, 'sys.get_eig')
            error_handling.positive_real(ms, 'ms')
            error_handling.positive_real(fs, 'fs')
            error_handling.lims(lims)
            fig, ax2 = plt.subplots()
            ax2 = plt.gca()
            if lims is None:
                ax2.set_ylim([-0.1, 1.1])
                ind = np.ones(self.sys.lat.sites, bool)
            else:
                ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
                ax2.set_ylim([lims[0]-0.1, lims[1]+0.1])
        else:
            ax2 = plt.twinx()
        error_handling.empty_ndarray(self.sys.pola, 'sys.get_pola')
        error_handling.tag(tag_pola, self.sys.lat.tags)
        x = np.arange(self.sys.lat.sites)
        i_tag = self.sys.lat.tags == tag_pola
        ax2.plot(x[ind], np.ravel(self.sys.pola[ind, i_tag]), 'or', markersize=(4*ms)//5)
        ylabel = '$<' + tag_pola.upper() + '|' + tag_pola.upper() + '>$'
        ax2.set_ylabel(ylabel, fontsize=fs, color='red')
        ax2.set_ylim([-0.1, 1.1])
        ax2.set_xlim(-0.5, x[ind][-1]+0.5)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
        for label in ax2.get_yticklabels():
            label.set_color('r')
        return fig, ax2

    def ipr(
        self,
        fig: Figure | None = None,
        ax1: Axes | None = None,
        ms: float = 10,
        fs: float = 20,
        lims: tuple[float, float] | None = None,
        ind: NDArray | None = None,
    ) -> tuple[Figure, Axes]:
        '''
        Plot Inverse Participation Ration.

        :param fig: Figure. Default value None. (used by the method spectrum).
        :param ax1: Axis. Default value None. (used by the method spectrum).
        :param ms: Positive Float. Default value 10. Markersize.
        :param fs: Positive Float. Default value 20. Fontsize.
        :param lims: List. lims[0] energy min, lims[1] energy max.
        :param ind: List. Default value None. List of indices. (used in the method spectrum).

        :returns:
            * **fig** -- Figure.
        '''
        if fig is None:
            error_handling.sys(self.sys)
            error_handling.empty_ndarray(self.sys.ipr, 'sys.get_ipr')
            error_handling.positive_real(ms, 'ms')
            error_handling.positive_real(fs, 'fs')
            error_handling.lims(lims)
            fig, ax2 = plt.subplots()
            ax2 = plt.gca()
            if lims is None:
                ind = np.ones(self.sys.lat.sites, bool)
            else:
                ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
        else:
            ax2 = plt.twinx()
        error_handling.empty_ndarray(self.sys.ipr, 'sys.get_ipr')
        x = np.arange(self.sys.lat.sites)
        ax2.plot(x[ind], self.sys.ipr[ind], 'or', markersize=(4*ms)//5)
        ax2.set_ylabel( 'IPR' , fontsize=fs, color='red')
        ax2.set_xlim(-0.5, x[ind][-1]+0.5)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
        for label in ax2.get_yticklabels():
            label.set_color('r')
        return fig, ax2

    def petermann(
        self,
        fig: Figure | None = None,
        ax1: Axes | None = None,
        ms: float = 10,
        fs: float = 20,
        lims: tuple[float, float] | None = None,
        ind: NDArray | None = None,
    ) -> tuple[Figure, Axes]:
        '''
        Plot Peterman factor.

        :param fig: Figure. Default value None. (used by the method spectrum).
        :param ax1: Axis. Default value None. (used by the method spectrum).
        :param ms: Positive Float. Default value 10. Markersize.
        :param fs: Positive Float. Default value 20. Fontsize.
        :param lims: List. lims[0] energy min, lims[1] energy max.
        :param ind: List. Default value None. List of indices. (used in the method spectrum).

        :returns:
            * **fig** -- Figure.
        '''
        if fig is None:
            error_handling.sys(self.sys)
            error_handling.empty_ndarray(self.sys.ipr, 'sys.get_petermann')
            error_handling.positive_real(ms, 'ms')
            error_handling.positive_real(fs, 'fs')
            error_handling.lims(lims)
            fig, ax2 = plt.subplots()
            ax2 = plt.gca()
            if lims is None:
                ind = np.ones(self.sys.lat.sites, bool)
            else:
                ind = (self.sys.en > lims[0]) & (self.sys.en < lims[1])
        else:
            ax2 = plt.twinx()
        error_handling.empty_ndarray(self.sys.ipr, 'sys.get_ipr')
        x = np.arange(self.sys.lat.sites)
        ax2.plot(x[ind], self.sys.petermann[ind], 'or', markersize=(4*ms)//5)
        ax2.set_ylabel( 'K' , fontsize=fs, color='red')
        ax2.set_xlim(-0.5, x[ind][-1]+0.5)
        for tick in ax2.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fs)
        for label in ax2.get_yticklabels():
            label.set_color('r')
        return fig, ax2

    def spectrum_complex(self, ms: float = 10., fs: float = 20., lims: tuple[float, float] | None = None) -> Figure:
        '''
        Plot complex value eigenenergies, real part (blue circles),
        and imaginary part (red circles).

        :param ms: Positive Float. Default value 20. Markersize.
        :param fs: Positive Float. Default value 20. Font size.
        :param lims: List. lims[0] energy min, lims[1] energy max.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.en, 'sys.get_eig')
        error_handling.positive_real(ms, 'ms')
        error_handling.positive_real(fs, 'fs')
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
        ax1.set_ylabel('Re '+r'$E_n$'+',    Im '+r'$E_n$', fontsize=fs)
        for label in ax1.xaxis.get_majorticklabels():
            label.set_fontsize(fs)
        for label in ax1.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        xa = ax1.get_xaxis()
        ax1.set_xlim([x[ind][0]-0.1, x[ind][-1]+0.1])
        xa.set_major_locator(plt.MaxNLocator(integer=True))
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def intensity_1d(
        self, intensity: NDArray, ms: float = 20., lw: float = 2., fs: float = 20.,
        title: str = r'$|\psi^{(j)}|^2$',
    ) -> Figure:
        '''
        Plot intensity for 1D lattices.

        :param intensity: np.array. Field intensity.
        :param ms: Positive Float. Default value 20. Markersize.
        :param lw: Positive Float. Default value 2. Linewith, connect sublattice sites.
        :param fs: Positive Float. Default value 20. Font size.
        :param title: String. Default value 'Intensity'. Figure title.
        '''
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.empty_ndarray(self.sys.lat.coor, 'sys.get_lattice')
        error_handling.positive_real(ms, 'ms')
        error_handling.positive_real(lw, 'lw')
        error_handling.positive_real(fs, 'fs')
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
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def intensity_disk(
        self,
        intensity: NDArray,
        s: float = 200.,
        fs: float = 20.,
        lims: tuple[float, float] | None = None,
        figsize: tuple[float, float] | None = None,
        title: str = r'$|\psi|^2$',
    ) -> Figure:
        r'''
        Plot the intensity. Colormap with identical disk shape.

        :param intensity: np.array.Field intensity.
        :param s: Default value 200. Disk size.
        :param fs: Default value 20. Font size.
        :param lims: List. Default value None. Colormap limits.
        :param figsize: Tuple. Default value None. Figure size.
        :param title: String. Default value '$|\psi_n|^2$'. Title.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.lat.coor, 'sys.get_lattice')
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.positive_real(s, 's')
        error_handling.positive_real(fs, 'fs')
        error_handling.tuple_2elem(figsize, 'figsize')
        error_handling.string(title, 'title')
        fig, ax = plt.subplots(figsize=figsize)
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
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def intensity_area(
        self,
        intensity: NDArray,
        s: float = 1000.,
        lw: float = 1.,
        fs: float = 20.,
        plt_hop: bool = False,
        figsize: tuple[float, float] | None = None,
        title: str = r'$|\psi|^2$',
    ) -> Figure:
        r'''
        Plot the intensity. Intensity propotional to disk shape.

        :param intensity: np.array. Intensity.
        :param s: Positive Float. Default value 1000.
            Circle size given by s * intensity.
        :param lw: Positive Float. Default value 1. Hopping linewidths.
        :param fs: Positive Float. Default value 20. Fontsize.
        :param plt_hop: Boolean. Default value False. Plot hoppings.
        :param figsize: Tuple. Default value None. Figure size.
        :param title: String. Default value '$|\psi_{ij}|^2$'. Figure title.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.sys.lat.coor, 'sys.get_lattice')
        error_handling.ndarray(intensity, 'intensity', self.sys.lat.sites)
        error_handling.positive_real(s, 's')
        error_handling.positive_real(fs, 'fs')
        error_handling.boolean(plt_hop, 'plt_hop')
        error_handling.tuple_2elem(figsize, 'figsize')
        error_handling.string(title, 'title')
        fig, ax = plt.subplots()
        ax.set_xlabel('$i$', fontsize=fs)
        ax.set_ylabel('$j$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        if plt_hop:
            plt.plot([self.sys.lat.coor['x'][self.sys.hop['i'][:]], 
                             self.sys.lat.coor['x'][self.sys.hop['j'][:]]],
                            [self.sys.lat.coor['y'][self.sys.hop['i'][:]],
                             self.sys.lat.coor['y'][self.sys.hop['j'][:]]],
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
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def butterfly(
        self,
        betas: NDArray,
        butterfly: NDArray,
        lw: float = 1.,
        fs: float = 20.,
        lims: tuple[float, float] | None = None,
        title: str = '',
    ) -> Figure:
        '''
        Plot energies depending on a parameter.

        :param betas: np.array. Parameter values.
        :param butterfly: np.array. Eigenvalues.
        :param lw: Positive Float. Default value 1. Hopping linewidths.
        :param fs: Positive Float. Default value 20. Fontsize.
        :param lims: List, lims[0] energy min, lims[1] energy max.
        :param title: Default value ''. Figure title.
        '''
        error_handling.ndarray_empty(betas, 'betas')
        error_handling.ndarray_empty(butterfly, 'butterfly')
        error_handling.positive_real(lw, 'lw')
        error_handling.positive_real(fs, 'fs')
        error_handling.lims(lims)
        error_handling.string(title, 'title')
        i_beta_min = np.argmin(np.abs(betas))
        if lims is None:
            lims = [butterfly[i_beta_min, 0], butterfly[i_beta_min, -1]]
        ind_en = np.argwhere((butterfly[i_beta_min, :] > lims[0]) & 
                                            (butterfly[i_beta_min, :] < lims[1]))
        ind_en = np.ravel(ind_en)
        fig, ax = plt.subplots()
        plt.title('Energies depending on strain', fontsize=fs)
        plt.xlabel(r'$\beta/\beta_{max}$', fontsize=fs)
        plt.ylabel('$E$', fontsize=fs)
        ax.set_title(title, fontsize=fs)
        plt.yticks(np.arange(lims[0], lims[1]+1, (lims[1]-lims[0])/4), fontsize=fs)
        plt.ylim(lims)
        beta_max = max(self.sys.betas)
        plt.xticks([-beta_max, -0.5*beta_max, 0, 
                        0.5*beta_max, beta_max], fontsize=fs)
        ax.set_xticklabels(('-1', '-1/2', '0', '1/2', '1'))
        plt.xlim([betas[0], betas[-1]])
        for i in ind_en:
            plt.plot(betas, butterfly[:, i], 'b', lw=lw)
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def show(self) -> None:
        """
        Emulate Matplotlib method plt.show().
        """
        plt.show()


# Backward-compatible lowercase alias (pre-0.2 API).
plot = Plot
