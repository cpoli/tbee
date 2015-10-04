'''
    **liebTB** package models a line centered squared lattice (Lieb lattice)
    with real and/or complex valued hoppings and onsite energies.
    Ponctual and dimerization sefects can be implemented.
    
    **liebTB** is available at https://github.com/cpoli/liebTB
'''

import numpy as np
import sys, os
from scipy.optimize import fsolve
import numpy.random as rand
from math import cos, sin, sqrt, pi
from cmath import exp as cexp  
from mpl_toolkits.mplot3d import axes3d
from matplotlib.legend_handler import HandlerLine2D
from latticeTB import *
from plotTB import *
from eigTB import *
from propagationTB import *

class eigLieb(eigTB):
    '''
    Child of the class **eigTB**. Dedicated to the Lieb lattice. 
    '''
    def __init__(self, coor, tags):
        eigTB.__init__(self, coor, tags)
        self.alpha = 0  # disorder strength
        self.c = 0  # next nearest hoppings strength 

    def set_hop_uni(self, t, c=0.):
        '''
        Get the nearest hoppings of the Lieb lattice and set their values to :math:`t`
        and, get the next nearest hoppings (diagonal hoppings)
        and set their values to :math:`c t`.

        :param t: Nearest hopping value.
        :param c: Next nearest hopping strength.
        '''
        self.c = c
        dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
        dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
        dis = dif_x **2 + dif_y **2
        ind = np.argwhere(dis == 1.)
        self.hop = np.zeros(len(ind[:, 0]), {'names': ['i', 'j', 't'], 
                                                              'formats': ['u4', 'u4', 'c16']})
        self.hop['i'] = ind[:, 0]
        self.hop['j'] = ind[:, 1]
        self.hop['t'] = 0.5 * t
        ind_nn = np.argwhere(dis == 2)
        hop_nn = np.zeros(len(ind_nn[:, 0]), {'names':['i', 'j', 't'],
                                                                   'formats':['u4', 'u4', 'c16']})
        hop_nn['i'] = ind_nn[:, 0]
        hop_nn['j'] = ind_nn[:, 1]
        hop_nn['t'] = 0.5 * c * t
        self.hop = np.concatenate([self.hop, hop_nn])

    def set_params(self):
        '''
        Set parameters used to store figures.

        :param on:  Array. Sublattices onsite energies.
        '''
        self.param = {'dis' :self.c, 'dis' :self.apha}
        print(self.param)

    def set_hop_alt(self, ta, tb, tc, td, c=0.):
        '''
        Get the nearest hoppings of the Lieb lattice with alternating hoppings
        and, the next nearest hoppings. As the  next nearest hoppings are 
        diagonal hoppings, the value of the diagonal between :math:`t_1` 
        and :math:`t_2` is given by :math:`c\sqrt{t_1^2+t_2^2}`.

        :param ta: Value of the math:`t_{ab}` hoppings.
        :param tb: Value of the math:`t_{ba}` hoppings.
        :param tc: Value of the math:`t_{ac}` hoppings.
        :param td: Value of the math:`t_{ca}` hoppings.
        :param c: Default value 0, strength of the diagonal hoppings.
        '''
        dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
        dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
        ind_x = np.argwhere((dif_x==1) & (dif_y==0))
        hop_x = np.zeros(len(ind_x[:, 0]), {'names': ['i', 'j', 't', 'tag'], 
                                                              'formats': ['u4', 'u4', 'c16', 'S2']})
        hop_x['i'] = ind_x[:, 0]
        hop_x['j'] = ind_x[:, 1]
        hop_x['tag'] = np.where(self.coor[ind_x[:, 0]]['tag'] == b'a', 'ab', 'ba')
        hop_x['t'] = np.where(hop_x['tag'] == b'ab', ta, tb)
        ind_y = np.argwhere((dif_x==0) & (dif_y==1))
        hop_y = np.zeros(len(ind_y[:, 0]), {'names':['i', 'j', 't', 'tag'], 
                                                              'formats':['u4', 'u4', 'c16', 'S2']})
        hop_y['i'] = ind_y[:, 0]
        hop_y['j'] = ind_y[:, 1]
        hop_y['tag'] = np.where(self.coor[ind_y[:, 0]]['tag'] == b'a', 'ac', 'ca')
        hop_y['t'] = np.where(hop_y['tag'] == b'ac', tc, td)
        self.hop = np.concatenate([hop_x, hop_y])
        #if c != 0:
        #    hop_nn = self.set_hop_alt_nn(c)
         #   self.hop = np.concatenate([self.hop, hop_nn])

    def set_hop_alt_nn(self, c):
        '''
        Get the next nearest hoppings if the Lieb lattice has alternating hoppings. 
        '''
        self.c = c
        dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
        dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
        ind = np.argwhere((dif_x == -1) & (dif_y == 1))
        n = len(ind[:, 0])
        ind_b = ind[self.coor['tag'][ind[:, 0]] == b'b']
        ind_c = ind[self.coor['tag'][ind[:, 0]] == b'c']
        print('ind_b', ind_b)
        print('ind_c', ind_c)
        hop_nn = np.zeros(2*n, {'names': ['i', 'j', 't', 'tag'],
                                      'formats': ['u4', 'u4', 'c16', 'S2']})
#        hop_pp = np.zeros(n, {'names': ['i', 'j', 't', 'tag'],
#                                     'formats': ['u4', 'u4', 'c16', 'S2']})
        hop_nn['i'][:n] = ind[:, 0]
        hop_nn['j'][:n] = ind[:, 1]
        for i, j in zip(ind_b[:,0], ind_b[:,1]):
            td = np.sqrt(self.hop['t'][(self.hop['i'] == i-1) & (self.hop['j'] == i)] ** 2 +
		     		  self.hop['t'][(self.hop['i'] == i-1) & (self.hop['j'] == j)] ** 2)
            hop_nn['t'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = c*td
            hop_nn['tag'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = b'-+'
        for i, j in zip(ind_c[:,0], ind_c[:,1]):
            td = np.sqrt(self.hop['t'][(self.hop['i'] == i) & (self.hop['j'] == j+1)] ** 2 +
                      self.hop['t'][(self.hop['i'] == j) & (self.hop['j'] == j+1)] ** 2)
            hop_nn['t'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = c*td
            hop_nn['tag'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = b'+-'
        ind = np.argwhere((dif_x == 1) & (dif_y == 1))
        n = len(ind[:, 0])
        ind_b = ind[self.coor['tag'][ind[:, 0]] == b'b']
        ind_c = ind[self.coor['tag'][ind[:, 0]] == b'c']
        ind_ord = ind[self.coor['tag'][ind[:, 0]] == b'b']
        hop_nn['i'][n:] = ind[:, 0]
        hop_nn['j'][n:] = ind[:, 1]
        for i, j in zip(ind_b[:,0], ind_b[:,1]):
            td = np.sqrt(self.hop['t'][(self.hop['i'] == i) & (self.hop['j'] == i+1)] ** 2 +
                      self.hop['t'][(self.hop['i'] == i+1) & (self.hop['j'] == j)] ** 2)
            hop_nn['t'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = c * td
            hop_nn['tag'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = b'++'
        for i, j in zip(ind_c[:,0], ind_c[:,1]):
            td = np.sqrt(self.hop['t'][(self.hop['i'] == j-1) & (self.hop['j'] == j)] ** 2 +
                      self.hop['t'][(self.hop['i'] == i) & (self.hop['j'] == j-1)] ** 2)
            hop_nn['t'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = c * td
            hop_nn['tag'][(hop_nn['i'] == i) & (hop_nn['j'] == j)]  = b'--'
        self.hop = np.concatenate([self.hop, hop_nn])

    def set_dim_defect(self, dim_x=-1, dim_y=-1):
        '''
        Set dimerization defects along :math:`x` and/or :math:`y`.

        :param dim_x: Default value -1 (no defect). Dimerization along :math:`x` site index.
        :param dim_x: Default value -1 (no defect). Dimerization along :math:`y` site index.
        '''
        ta = self.hop['t'][self.hop['tag'] == b'ab'][0]
        tb = self.hop['t'][self.hop['tag'] == b'ba'][0]
        tc = self.hop['t'][self.hop['tag'] == b'ac'][0]
        td = self.hop['t'][self.hop['tag'] == b'ca'][0]
        if dim_x > 1:
            self.hop['t'][(self.hop['tag'] == b'ab') & 
		      (self.coor['x'][self.hop['i']] >= dim_x)] = tb
            self.hop['t'][(self.hop['tag'] == b'ba') & 
		      (self.coor['x'][self.hop['i']]  >= dim_x)] = ta
        if dim_y > 1:
            self.hop['t'][(self.hop['tag'] == b'ac') & 
		      (self.coor['y'][self.hop['i']] >= dim_y)] = td
            self.hop['t'][(self.hop['tag'] == b'ca') &
		      (self.coor['y'][self.hop['i']]  >= dim_y)] = tc

    def disorder_uniform(self, nx, ny, alpha=0):
        '''
        Set a disorder uniform along math:`y` for math:`t_{ab}`  and math:`t_{ba}`,  
        and uniform along math:`x` for math:`t_{ac}`  and math:`t_{ca}`.

        :param alpha: Default value 0. Stength of the disorder.

        .. note ::

	This disorder preserves the zero mode.
        '''
        self.alpha = alpha
        rand_x = (1.+alpha*rand.uniform(-1, 1., nx-1))
        rand_y = (1.+alpha*rand.uniform(-1, 1., ny-1))
        for y in range(0, ny, 2):
            self.hop['t'][((self.hop['tag'] == b'ab') | (self.hop['tag'] == b'ba')) &
		       (self.coor['y'][self.hop['i']] == y)] *= rand_x
        for x in range(0, nx, 2):
            self.hop['t'][((self.hop['tag'] == b'ac') | (self.hop['tag'] == b'ca')) &
		       (self.coor['x'][self.hop['i']] == x)] *= rand_y

    def disorder_pair(self, nx, ny, alpha=0):
        '''
        Set a non uniform disorder.

        :param alpha: Default value 0. Stength of the disorder.

        .. note ::

	This disorder preserves the zero mode.
        '''
        self.alpha = alpha
        for y in range(0, ny, 2):
            rand_x = (1+alpha*rand.uniform(-1, 1., nx-1))
            rand_x[1::2] = rand_x[0::2]
            self.hop['t'][((self.hop['tag'] == b'ab') | (self.hop['tag'] == b'ba')) & 
		       (self.coor['y'][self.hop['i']] == y)] *= rand_x
        for x in range(0, nx, 2):
            rand_y = (1+alpha*rand.uniform(-1, 1., ny-1))
            rand_y[1::2] = rand_y[0::2]
            self.hop['t'][((self.hop['tag'] == b'ac') | (self.hop['tag'] == b'ca')) &
		       (self.coor['x'][self.hop['i']] == x)] *= rand_y

    def disorder_plaquette(self, nx, ny, alpha=0.):
        '''
        Set a non uniform disorder.

        :param alpha: Default value 0. Stength of the disorder.

        .. note ::

	This disorder preserves the zero mode.
        '''
        self.alpha = alpha
        lx, ly = 2, 2
        self.hop['t'] *= 1+ alpha * rand.uniform(-1, 1., len(self.hop['tag']))
        for i, j in [(i, j) for i in range((nx-1)//2) for j in range((ny-1)//2)]: 
            hop_p = self.hop[(self.coor['x'][self.hop['i']]>i*lx-0.1) &
			    (self.coor['x'][self.hop['i']]<(i+1)*lx+0.1) &
			    (self.coor['y'][self.hop['i']]>j*ly-0.1) &
			    (self.coor['y'][self.hop['i']]<(j+1)*ly+0.1) &
			    (self.coor['x'][self.hop['j']]>i*lx-0.1) &
		                (self.coor['x'][self.hop['j']]<(i+1)*lx+0.1) &
			    (self.coor['y'][self.hop['j']]>j*ly-0.1) &
			    (self.coor['y'][self.hop['j']]<(j+1)*ly+0.1)]
            self.hop['t'][(self.hop['i'] == hop_p['i'][7]) & (self.hop['j'] == hop_p['j'][7])] = \
		       hop_p['t'][0]*hop_p['t'][3]*hop_p['t'][5]*hop_p['t'][6] / \
		       (hop_p['t'][1]*hop_p['t'][2]*hop_p['t'][4])

    def magnetic_field(self, nx):
        '''
        Add Pierls phases to the hoppings.

        :param nx: Number of sites along math:`x`. (Input of the class **latticeTB**)

        .. note ::

            This disorder preserves the zero mode.
        '''
        for x in range(0, nx, 2):
            self.hop['t'][((self.hop['tag'] == b'ac') | (self.hop['tag'] == b'ca')) 
                            & (self.coor['x'][self.hop['i']] == x)] *= \
		          np.exp(1j * pi * flux * x//2)


class plotLieb(plotTB):
    '''
    Plot the output of the class **latticeTB** and the class **eigTB**.
    '''
    def __init__(self, coor, tags):
        plotTB.__init__(self, coor, tags)

    def plt_lieb_hop(self, hop, nx, ny, ms=20):
        '''
        Plot the coordinates of  the Lieb lattice sites in hoppings space.

        :param hop: Hoppings.
        :param nx: Number of sites along :math:`x`.
        :param ny: Number of sites along :math:`y`.	
        :param ms: Default value 20. Markersize.
        '''
        lx, ly = 2, 2
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ax.axis('off')
        coor_hop = np.copy(self.coor)
        # coordinates at the bottom 
        coor_hop['x'][1: nx] = np.cumsum(hop['t'][: nx-1].real)
        # coordinates at the left 
        ty = hop['t'][((hop['tag'] == b'ac') | (hop['tag'] == b'ca')) & \
		  ((self.coor['x'][hop['i']] == 0) & (self.coor['x'][hop['j']] == 0))].real
        coor_hop['y'][(self.coor['x'] == 0) & (self.coor['y'] > 0)] = np.cumsum(ty)
        # loop over each plaquette
        for i, j in [(i, j) for i in range(nx//2) for j in range(ny//2)]: 
            ind_r = np.where((self.coor['x']>i*lx-0.1) & (self.coor['x']<(i+1)*lx+0.1) &
                                        (self.coor['y']>j*ly-0.1) & (self.coor['y']<(j+1)*ly+0.1))
            ind_r = np.ravel(ind_r)
            cond_h = (self.coor['x'][hop['i']]>i*lx-0.1) & (self.coor['x'][hop['i']]<(i+1)*lx+0.1) & \
		    (self.coor['y'][hop['i']]>j*ly-0.1) & (self.coor['y'][hop['i']]<(j+1)*ly+0.1) & \
		    (self.coor['x'][hop['j']]>i*lx-0.1) & (self.coor['x'][hop['j']]<(i+1)*lx+0.1) & \
		    (self.coor['y'][hop['j']]>j*ly-0.1) & (self.coor['y'][hop['j']]<(j+1)*ly+0.1)
            hop_p = hop['t'][cond_h].real
            tx0, tx1, tx2, tx3 = hop_p[0], hop_p[1], hop_p[2], hop_p[3]
            ty0, ty2, ty1, ty3 = hop_p[4], hop_p[5], hop_p[6], hop_p[7] 

            p = (coor_hop['x'][ind_r[5]], coor_hop['y'][ind_r[5]],
                   coor_hop['x'][ind_r[2]], coor_hop['y'][ind_r[2]], tx2+tx3, ty2+ty3)
            if np.sqrt((p[0]-p[2])**2+(p[1]-p[3])**2) > p[4]+p[5]:
                print('\n\nWARNING WARNING WARNING')
                print('Lattice construction in hopping space is not possible.')
                raise Exception('\n\nRun again or decrease the strength '
				'of the disorder.\n\n')
            ang_x, ang_y = fsolve(equations, (0., 0.5*pi), (p,))
            coor_hop['x'][ind_r[4]] = coor_hop['x'][ind_r[2]] + ty2*cos(ang_y)
            coor_hop['y'][ind_r[4]] = coor_hop['y'][ind_r[2]] + ty2*sin(ang_y)
            coor_hop['x'][ind_r[6]] = coor_hop['x'][ind_r[5]] + tx2*cos(ang_x)
            coor_hop['y'][ind_r[6]] = coor_hop['y'][ind_r[5]] + tx2*sin(ang_x) 
            coor_hop['x'][ind_r[7]] = coor_hop['x'][ind_r[4]] + ty3*cos(ang_y) 
            coor_hop['y'][ind_r[7]] = coor_hop['y'][ind_r[4]] + ty3*sin(ang_y)
        # plot bonds
        for i in range(len(hop)): 
            plt.plot([coor_hop['x'][hop['i'][i]], coor_hop['x'][hop['j'][i]]],
                        [coor_hop['y'][hop['i'][i]], coor_hop['y'][hop['j'][i]]],
                 	'k', lw=1)
        # plot sites
        for t, c in zip(self.tags, self.colors):
            plt.plot(coor_hop['x'][coor_hop['tag'] == t],
		coor_hop['y'][coor_hop['tag'] == t],
		'o', color=c, ms=ms, markeredgecolor='none')
        xlim = [np.min(coor_hop['x'])-1.5, np.max(coor_hop['x'])+1.5]
        ylim = [np.min(coor_hop['y'])-1.5, np.max(coor_hop['y'])+1.5]
        ax.set_aspect('equal')
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.draw()
        return fig

    def plt_dispersion(self, ta, tb, tc, td, fs=20, c=0., N=100):
        '''
        Plot the relation dispersion of the Lieb lattice.

        :param ta: Hoppings.
        :param nx: Number of sites along :math:`x`.
        :param ny: Number of sites along :math:`y`.	
        :param ms: Default value 20. Markersize.
        '''
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        Ek = np.zeros((N, N, 3))
        kx = np.linspace(0, 2*pi, N)
        ky = np.linspace(0, 2*pi, N)
        if not c:
            kx , ky = np.meshgrid(kx , ky)
            f = ta + tb * np.exp(-1j * kx)
            g = tc + td * np.exp(-1j * ky)	
            en = np.sqrt(np.abs(f) **2 + np.abs(g) **2)
            Ek[:, :, 0] = - en
            Ek[:, :, 1] = 0
            Ek[:, :, 2] = en
        else:
            f = ta + tb * np.exp(-1j * kx)
            g = tc + td * np.exp(-1j * ky)
            t1 = c * sqrt(ta **2+tc **2)
            t2 = c * sqrt(tb **2+tc **2)
            t3 = c * sqrt(ta **2+td **2)
            t4 = c * sqrt(tb **2+td **2)
            for i, j in [(i, j) for i in range(N) for j in range(N)]:
                l = t1 + t2 * cexp(1j * kx[i]) + t3 * cexp(-1j * ky[j]) \
                          + t4 * cexp(1j *( kx[i] - ky[j]))
                h = np.array([[0., f[i], g[j]], 
                                      [f[i].conj(), 0., l],
                                      [g[j].conjugate(), l.conjugate(), 0.]])
                Ek[i, j, :] = LA.eigvalsh(h)
            kx , ky = np.meshgrid(kx , ky)
        alpha = 0.8
        ax.plot_surface(kx, ky, Ek[:, :, 0], linewidth=0.1, antialiased=True, 
		          shade=True, alpha=alpha, color='b')
        ax.plot_surface(kx, ky, Ek[:, :, 1], linewidth=0.1, antialiased=True, 
		          shade=True, alpha=alpha, color='g')
        ax.plot_surface(kx, ky, Ek[:, :, 2], linewidth=0.1, antialiased=True, 
			shade=True, alpha=alpha, color='r')
        plt.title('$E(\mathbf{k})$', fontsize=fs)
        ax.set_xlabel('$k_x$', fontsize=fs)
        ax.set_ylabel('$k_y$', fontsize=fs)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([0])
        plt.draw()
        return fig


def equations(p, param):
    '''
    Get the intersection of the points at the extremity of two segments.
    Used to build up the lattice in hopping space.
    '''
    (ta, tb) = p
    xa, ya, xb, yb, la, lb = param
    return (xa+la*cos(ta)-xb-lb*cos(tb),
            ya+la*sin(ta)-yb-lb*sin(tb))

"""
def basic_animation(frames=10, interval=30, s = 300):
    '''Plot a basic sine wave with oscillating amplitude'''
    fig = plt.figure()
    x = np.linspace(0, 10, 10)
    y = np.linspace(0, 10, 10)
    c = np.linspace(-5, 5, 10)
    
    ticks = [-np.max(c), np.max(c)]
    scat = plt.scatter(x, y, c=c, s=s, vmin=ticks[0], vmax=ticks[1],
                     cmap=plt.get_cmap('seismic'))
    frame = plt.gca()
    frame.axes.get_xaxis().set_ticks([])
    frame.axes.get_yaxis().set_ticks([])
    cbar = fig.colorbar(scat, ticks=ticks)
    cbar.ax.set_yticklabels(['', ''])
    def init():
        scat.set_array(c)
        return scat,

    def animate(i):
        scat.set_array(c)
        return scat,

    return animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frames, interval=interval)
"""
'''
# Lieb lattice inputs
nx, ny = 9, 9
ri = [[0, 0], [1, 0], [0, 1]]
pv = [[2, 2], [2, 2], [2, 2]]
tags = np.array([b'a', b'b', b'c'])
lat_lieb = latticeTB(tags=tags, ri=ri, pv=pv)
lat_lieb.get_lattice(nx=nx, ny=ny)
plt_lieb = plotTB(lat_lieb.coor, tags)
fig_lat = plt_lieb.plt_lattice(ms=20, fs=25, numero=True)
eig_lieb = eigLieb(lat_lieb.coor, tags)

# UNIFORM HOPPINGS
t = 1. # uniform hopping
on = [0j, -.5j, -0.5j] # onsite energies
eig_lieb.set_hop_uni(t=t)
eig_lieb.set_onsite(on=on)
eig_lieb.get_ham()
eig_lieb.get_eig()
fig_spec = plt_lieb.plt_spec(eig_lieb.en, eig_lieb.pola, pola_tag=b'a', ms=20)
zero_mode = eig_lieb.get_state_pola(b'a')
'''
'''
# Lieb lattice inputs
ta, tb = 2., 1.
tc, td = 2., 1.
nx, ny = 9, 9
ri = [[0, 0], [1, 0], [0, 1]]
pv = [[2, 2], [2, 2], [2, 2]]
t = 1
tags = np.array([b'a', b'b', b'c'])
lat_lieb = latticeTB(tags=tags, ri=ri, pv=pv)
lat_lieb.get_lattice(nx=nx, ny=ny)
plt_lieb = plotLieb(lat_lieb.coor, tags)
plt_lieb.plt_dispersion(ta=ta, tb=tb, tc=tc, td=td, N=200, fs=20, c=0.5)
'''
'''
fig_lat = plt_lieb.plt_lattice(ms=20, fs=25, numero=True)
eig_lieb = eigLieb(lat_lieb.coor, tags)
eig_lieb.set_hop_uni(t=t)
'''
nx, ny = 9, 9
ri = [[0, 0], [1, 0], [0, 1]]
pv = [[2, 2], [2, 2], [2, 2]]
tags = np.array([b'a', b'b', b'c'])
lat_lieb = latticeTB(tags=tags, ri=ri, pv=pv)
lat_lieb.get_lattice(nx=nx, ny=ny)
plt_lieb = plotTB(lat_lieb.coor, tags)
eig_lieb = eigLieb(lat_lieb.coor, tags)

ta, tb = 2., 1.
tc, td = 2., 1.
onsite = [0j, -1j, -1j]
eig_lieb.set_hop_alt(ta=ta, tb=tb, tc=tc, td=td, c=0.)
eig_lieb.set_dim_defect(dim_x=nx//2, dim_y=nx//2)
eig_lieb.set_hop_alt_nn(c=1)
#eig_lieb.disorder_plaquette(nx=nx, ny=ny, alpha=0.1)
eig_lieb.set_onsite(onsite)
eig_lieb.get_ham()
eig_lieb.get_eig()
#fig_spec = plt_lieb.plt_spec(eig_lieb.en, eig_lieb.pola, ms=20, pola_tag=b'a')
zero_mode = eig_lieb.get_state_pola(b'a')
#fig_state = plt_lieb.plt_intensity(zero_mode)
fig_lat = plt_lieb.plt_lattice(ms=20, fs=25, label=False, hop=eig_lieb.hop)
#fig_lat_hop = plt_lieb.plt_lieb_hop(hop=eig_lieb.hop, nx=nx, ny=ny, ms=15)
#save = saveFigTB(sys=eig_lieb, hop_tags=['ta', 'tb', 'tc', 'td'], hop=[ta, tb, tc, td], 
#					name_dir='lieb', ext='png')
#save.save_fig(fig_lat_hop, 'lattice_dimXY')
#save.save_fig(fig_spec, 'spec_dimXY')
#save.save_fig(fig_state, 'mode_dimXY')



'''
xy = np.loadtxt('lattice_real_dimDefx_4_dimDefy_4_dis_gen_x0,35_y0,35.txt')
x = dist_to_hop_mwa(xy[:, 0])
y = dist_to_hop_mwa(xy[:, 1])
lat_lieb.coor['x'] = x
lat_lieb.coor['y'] = y
print(lat_lieb.coor)
fig = plt.figure()
plt.plot(x, y, 'o')

#plt_lieb.plt_lattice(ms=20, fs=25, numero=False)
'''

plt.show()

