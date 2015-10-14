import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from JSAnimation import IPython_display
import os


class propTB():
    '''
    Get lattice time evolution. Time dependent Schrödinger equation solved by
    Crank-Nicolson method.
    '''
    def __init__(self, lat, steps, dz):
        '''
        Plot the coordinates of  the Lieb lattice sites in hoppings space.

        :param lat: **latticeTB** class instance.
        :param steps: Number of steps.
        :param dz: Step.
        '''
        self.lat = lat
        self.steps = steps
        self.dz = dz
        self.prop = np.zeros((self.lat.sites, self.steps), 'c16')

    def get_prop(self, ham, psi_init, norm=True):
        '''
        Get the time evolution.

        :param ham: Tight-Binding Hamilonian.
        :param psi_init: Initial state.
        :param norm: Default value True. Normalize the norm to 1 at each step.
        '''
        psi_init = np.array([psi_init])
        self.prop[:, 0] = psi_init
        diag = 1j*np.ones(self.lat.sites, 'c16')
        A = (sparse.diags(diag, 0) - 0.5 * self.dz * ham).toarray()
        B = (sparse.diags(diag, 0) + 0.5 * self.dz * ham).toarray()
        mat = (np.dot(LA.inv(A), B))
        for i in range(1, self.steps):
            self.prop[:, i] = np.dot(mat, self.prop[:, i-1])
            if norm:
                self.prop[:, i] /= np.sqrt((np.abs(self.prop[:, i])**2).sum())

    def get_pump(self, hams, psi_init, norm=True):
        '''
        Get the time evolution under adiabatic pumpings.

        :param hams: Tight-Binding Hamilonians.
        :param psi_init: Initial state.
        :param norm: Default value True. Normalize the norm to 1 at each step.
        '''
        ham = np.array([hams])
        psi_init = np.array([psi_init])
        no, = hams.shape
        self.prop[:, 0] = psi_init
        diag = 1j*np.ones(self.lat.sites, 'c16')
        delta = self.steps //(no+1)
        A = (sparse.diags(diag, 0) - 0.5 * self.dz * hams[0]).toarray()
        B = (sparse.diags(diag, 0) + 0.5 * self.dz * hams[0]).toarray()
        mat = (np.dot(LA.inv(A), B))
        # before pumping
        for i in range(1, delta):
           self.prop[:, i] = np.dot(mat, self.prop[:, i-1])
           if norm:
               self.prop[:, i] /= np.sqrt((np.abs(self.prop[:, i])**2).sum())
      # pumping
        c = np.linspace(0, 1, delta)
        for j in range(0, no-1):
            for i in range(0, delta):
                ham = (1-c[i])*hams[j]+c[i]*hams[j+1]
                A = (sparse.diags(diag, 0) - 0.5 * self.dz * ham).toarray()
                B = (sparse.diags(diag, 0) + 0.5 * self.dz * ham).toarray()
                mat = (np.dot(LA.inv(A), B))
                self.prop[:, (j+1)*delta+i] = np.dot(mat, self.prop[:,  (j+1)*delta+i-1])
                if norm:
                    self.prop[:,  (j+1)*delta+i] /= \
                        np.sqrt((np.abs(self.prop[:,  (j+1)*delta+i])**2).sum())    
      # after pumping
        j = no
        delta_f = self.steps - no*delta
        for i in range(0, delta_f):
            self.prop[:,  j*delta+i] = np.dot(mat, self.prop[:,  j*delta+i-1])
            if norm:
                self.prop[:,  j*delta+i] /= np.sqrt((np.abs(self.prop[:,  j*delta+i])**2).sum())

    def plt_prop1d(self, fs=20):
        '''
        Plot time evolution for 1D systems. 

        :param fs: Default value 20. Fontsize.
        '''
        if not self.prop.any():
            raise Exception('\n\nRun method get_prop() or get_pump() first.\n')
        fig, ax = plt.subplots(figsize=(8, 6))
        prop = np.abs(self.prop[::-1, :]) **2
        color = self.prop_smooth1d(prop)
        plt.title('$|\psi(z)|^2$', fontsize=fs)
        plt.ylabel('n', fontsize=fs)
        plt.xlabel('z', fontsize=fs)
        vmin, vmax = 0., np.max(color[:, 0: self.lat.sites])
        
        extent = (0, self.steps*self.dz, 
                       -self.lat.sites//2+0.5, self.lat.sites-self.lat.sites//2+0.5)
        aspect = 'auto'
        interpolation = 'nearest'
        im = plt.imshow(color, cmap=plt.cm.hot, aspect=aspect,
                                  interpolation=interpolation, extent=extent,
                                  vmin=vmin, vmax=vmax)
        
        cbar = plt.colorbar(im, ticks=[vmin, vmax])
        cbar.ax.set_yticklabels(['0', 'max'], fontsize=fs)
        ya = ax.get_yaxis()
        ya.set_major_locator(plt.MaxNLocator(integer=True))
        return fig

    def prop_smooth1d(self, prop, a=14., no=40):
        r'''
        Smooth propagation for 1D systems.
        Perform Gaussian interpolation :math:`e^{-a(x-x_i)^2}`,

        :param prop: Propagation.
        :param a: Default value 10. Gaussian Parameter.
        :param no: Default value 40. Number of points of each Gaussians.

        :returns:
           * **smooth** -- Smoothed propagation.
        '''
        x = np.linspace(-0.5, 0.5, no)
        smooth = np.empty((self.lat.sites * no, self.steps))
        for iz in range(0, self.steps):
            for i in range(self.lat.sites):
                smooth[i*no: (i+1)*no, iz] = prop[i, iz] * np.exp(-a * x ** 2)
        return smooth

    def get_ani(self, s=300, fs=20):
        '''
        Get time evolution animation.

        :param s: Default value 300. Circle shape.
        :param fs: Default value 20. Fontsize.

        :returns:
          * **ani** -- Animation.
        '''
        if not self.prop.any():
            raise Exception('\n\nRun method get_prop() or get_pump() first.\n')
        if os.name == 'posix':
            blit = False
        else:
            blit = True
        color = self.prop.real
        fig = plt.figure()
        plt.xlim([self.lat.coor['x'][0]-.5, self.lat.coor['x'][-1]+.5])
        plt.ylim([self.lat.coor['y'][0]-.5, self.lat.coor['y'][-1]+.5])
        ticks = [-np.max(color), np.max(color)]
        scat = plt.scatter(self.lat.coor['x'], self.lat.coor['y'], c=color[:, 0],
                                   s=s, vmin=ticks[0], vmax=ticks[1],
                                   cmap=plt.get_cmap('seismic'))
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        frame.axes.get_yaxis().set_ticks([])
        cbar = fig.colorbar(scat, ticks=[ticks[0], 0, ticks[1]])
        cbar.ax.set_yticklabels(['', '0',''])
        def update(i, color, scat):
            scat.set_array(color[:, i])
            return scat,
        ani = animation.FuncAnimation(fig, update, frames=self.steps,
                                                  fargs=(color, scat), blit=blit, repeat=False)
        return ani

    def get_ani_nb(self, s=300, fs=20):
        '''
        Get time evolution animation for iPython notebooks.

        :param s: Default value 300. Circle shape.
        :param fs: Default value 20. Fontsize.

        :returns:
           * **ani** -- Animation.
        '''
        if not self.prop.any():
            raise Exception('\n\nRun method get_prop() or get_pump() first.\n')
        fig = plt.figure()
        color = self.prop.real
        fig = plt.figure()
        ticks = [-np.max(color), np.max(color)]
        scat = plt.scatter(self.lat.coor['x'], self.lat.coor['y'], c=color[:, 0],
                                    s=s, vmin=ticks[0], vmax=ticks[1],
                                   cmap=plt.get_cmap('seismic'))
        frame = plt.gca()
        frame.axes.get_xaxis().set_ticks([])
        frame.axes.get_yaxis().set_ticks([])
        cbar = fig.colorbar(scat, ticks=[ticks[0], 0, ticks[1]])
        cbar.ax.set_yticklabels(['', '0',''])
        def init():
             scat.set_array(color[:, 0])
             return scat,
        def update(i):
             scat.set_array(color[:, i])
             return scat,
        return animation.FuncAnimation(fig, update, init_func=init,
                                                           frames=self.steps, interval=120)

    def plt_prop_dimer(self, lw=5, fs=20):
        '''
        Plot time evolution for dimers.
        
        :param lw: Default value 5. Linewidth.
        :param fs: Default value 20. Fontsize.

        :returns:
           * **fig** -- Figure.
        '''
        if not self.prop.any():
            raise Exception('\n\nRun method get_prop() or get_pump() first.\n')
        color = ['b', 'r']
        fig, ax = plt.subplots()
        z = self.dz * np.arange(self.steps)
        for i, c in zip([0, 1], color):
            plt.plot(z, np.abs(self.prop[i, :])**2, c, lw=lw)
        plt.title('Intensity', fontsize=fs)
        plt.xlabel('$z$', fontsize=fs)
        plt.ylabel('$|\psi_j|^2$', fontsize=fs)
        plt.xlim([0, z[-1]])
        return fig


