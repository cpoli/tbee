import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import sys
import numpy.core.defchararray as npc
from math import sin, cos
PI = np.pi
import error_handling


def empty_array(arr):
    length = len(arr)
    if length:
        return np.delete(arr, range(length))
    else:
        return arr


class system():
    def __init__(self, lat):
        '''
        Solve the Tight-Binding eigenvalue problem of a lattice defined 
        by the class **lattice**.

        :param lat: **lattice** class instance.
        '''
        error_handling.lat(lat)
        self.lat = lat
        self.coor_hop = np.array([], dtype=[ ('x','f16'), ('y','f16'), ('tag','S1')])
        self.vec_hop = np.array([], dtype=[('d', 'f16'),  ('a', 'f16')]) # Hopping lengths and distances
        self.dist_uni = np.array([], 'f8')  # Different hopping lengths
        self.ang_uni = np.array([], 'f8')  # Different hopping angles
        self.hop = np.array([], dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                       ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')]) #  Hoppings
        self.onsite = np.array([], 'c16') #  Onsite energies
        self.ham = sparse.csr_matrix(([],([],[])), shape=(self.lat.sites, self.lat.sites))  # Hamiltonian
        self.en = np.array([], 'c16')  # Eigenenergies
        self.vn = np.array([], 'c16')  # Eigenvectors
        self.intensity = np.array([], 'f8')  # Intensities (|vn|**2)
        self.ipr = np.array([], 'f8')  # Inverse Participation Ratio (|vn|**4)
        self.pola = np.array([], 'f8')  # sublattices polarisation (|vn^{(S)}|**2)
        self.alpha = 0.  # hopping disorder strength
        self.alpha_onsite = 0.  # onsite disorder strength
        self.params = {}

    def clean_hopping(self):
        '''
        Clean structured array *hop*
        '''
        self.hop = np.array([], dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                       ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])

    def get_hopping(self):
        '''
        Get lengths and angles of the hoppings.
        '''
        error_handling.sites(self.lat.sites)
        dif_x = self.lat.coor['x'] - self.lat.coor['x'].reshape(self.lat.sites, 1)
        dif_y = self.lat.coor['y'] - self.lat.coor['y'].reshape(self.lat.sites, 1)
        dist = np.sqrt(dif_x ** 2 + dif_y ** 2).round(3)
        ang = (180 / PI * np.arctan2(dif_y, dif_x)).round(3)
        self.vec_hop = np.zeros(dist.shape, dtype=[('d', 'f16'),  ('a', 'f16')])
        self.vec_hop['d'] = dist
        self.vec_hop['a'] = ang
        self.dist_uni = np.unique(self.vec_hop['d'])
        self.ang_uni = np.unique(self.vec_hop['a'])

    def print_hopping(self, n=5):
        '''
        Print the distances and the angles of all hoppings.

        :param n: Positive integer. Print the first nth hopping
          distances and associated positive angles.
        '''
        error_handling.sites(self.lat.sites)
        self.get_hopping()
        nmax = len(self.dist_uni) - 1
        error_handling.positive_int_lim(n, 'n', nmax)
        print('\n{} different distances between sites:'.format(nmax))
        print('\nDistances between sites:')
        for i, d in enumerate(self.dist_uni[1:n+1]):
            if i == 0:
                hop_name = 'st'
            elif i == 1:
                hop_name = 'nd'
            elif i == 2:
                hop_name = 'rd'
            else:
                hop_name = 'th'
            print('{}{} hopping, length: {:.3f}'.format(i+1, hop_name, d))
            print('\twith positive angles:')
            positive_ang = self.vec_hop['a'][(self.vec_hop['d'] == d) &
                                                               (self.vec_hop['a'] >= 0.) &
                                                               (self.vec_hop['a'] < 180.)]
            print('\t', np.unique(positive_ang))

    def set_onsite(self, dict_onsite):
        '''
        Set onsite energies.

        :param on:  Array. Sublattice onsite energies.
        '''
        error_handling.sites(self.lat.sites)
        error_handling.set_onsite(dict_onsite, self.lat.tags)
        self.onsite = np.zeros(self.lat.sites, 'c16')
        for tag, on in dict_onsite.items():
            self.onsite[self.lat.coor['tag'] ==tag] = on

    def set_hopping(self, list_hop):
        '''
        Set lattice hoppings.

        :param list_hop: List of dictionaries, Dictionary with key a tuple:(n, 'ang') nth hopping,
          associated positive angle, and hopping value {val}.
        '''        
        error_handling.sites(self.lat.sites)
        self.get_hopping()
        error_handling.set_hopping(list_hop, len(self.dist_uni) - 1)
        list_n = np.unique([dic['n'] for dic in list_hop])
        for n in list_n:
            self.hop = np.delete(self.hop, np.where(self.hop['n'] == n))
        for n in list_n:
            ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[n] - 1e-3) &
                                       (self.vec_hop['d'] < self.dist_uni[n] + 1e-3))
            ind_up = ind[ind[:, 1] > ind[:, 0]]
            hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                               ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
            hop['n'] = n
            hop['i'] = ind_up[:, 0]
            hop['j'] = ind_up[:, 1]
            hop['ang'] = self.vec_hop['a'][ind_up[:, 0], ind_up[:, 1]]
            hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], self.lat.coor['tag'][ind_up[:, 1]])
            for dic in list_hop:
                if dic['n'] != n:
                    continue
                if len(dic) == 2:
                    hop['t'] = dic['t']
                elif len(dic) == 3 and 'ang' in dic:
                    hop['t'][hop['ang'] == dic['ang']] = dic['t']
                elif len(dic) == 3 and 'tag' in dic:
                    hop['t'][hop['tag'] == dic['tag']] = dic['t']
                else:
                    hop['t'][(hop['tag'] == dic['tag']) & (hop['ang'] == dic['ang'])] = dic['t']
            self.hop = np.concatenate([self.hop, hop])

    def  set_hopping_dis(self, alpha):
        '''
        Set uniform hopping disorder. 

        :param alpha: Number. Stength of the disorder.
        '''
        error_handling.empty_hop(self.hop)
        error_handling.number(alpha, 'alpha')
        self.hop['t'] *= 1. + alpha * rand.uniform(-1., 1., len(self.hop))
        self.alpha = alpha


    def set_onsite_dis(self, alpha):
        '''
        Set uniform onsite disorder. 

        :param alpha: Number. Stength of the disorder.
        '''
        error_handling.empty_onsite(self.onsite)
        error_handling.number(alpha, 'alpha')
        self.onsite *= 1. + alpha * rand.uniform(-1., 1., self.lat.sites)
        self.alpha_onsite = alpha

    def set_onsite_def(self, onsite_def):
        '''
        Set specific onsite energies.

        :param dict_ons_def:  Dictionary. key: site indices, val: onsite values. 
        '''
        error_handling.sites(self.lat.sites)
        error_handling.empty_onsite(self.onsite)
        error_handling.set_onsite_def(onsite_def, self.lat.sites)
        for i, o in onsite_def.items():
            self.onsite[i] = o

    def set_hopping_def(self, hopping_def):
        '''
        Set specific hoppings. 

        :param dict_hop_def:  Dictionary. key: hopping indices, val: hopping values. 
        '''
        error_handling.empty_hop(self.hop)
        error_handling.set_hopping_def(self.hop, hopping_def, self.lat.sites)
        for key, val in hopping_def.items():
            cond = (self.hop['i'] == key[0]) & (self.hop['j'] == key[1])
            self.hop['t'][cond] = val
            cond = (self.hop['j'] == key[0]) & (self.hop['i'] == key[1])
            self.hop['t'][cond] = val

    def change_hopping(self, list_hop, x_bottom_left=0, y_bottom_left=0):
        '''
        Change hopping values.

        :param dict_hop: Dictionary. key a tuple:(n, 'ang'} nth hopping,
          associated positive angle, and hopping value {val}.
        :param x_bottom_left: Real number. lower bound along:math:`x` 
        :param y_bottom_left: Real number. lower bound along:math:`y` 
        '''
        error_handling.empty_hop(self.hop)
        error_handling.set_hopping(list_hop, len(self.dist_uni) - 1)
        error_handling.real_number(x_bottom_left, 'x_bottom_left')
        error_handling.real_number(y_bottom_left, 'y_bottom_left')
        ind = (self.lat.coor['x'][self.hop['i']] >= x_bottom_left) & \
                 (self.lat.coor['y'][self.hop['i']] >= y_bottom_left) & \
                 (self.lat.coor['x'][self.hop['j']] >= x_bottom_left) & \
                 (self.lat.coor['y'][self.hop['j']] >= y_bottom_left)
        for dic in list_hop:
            for dic in list_hop:
                if len(dic) == 2:
                    self.hop['t'][ind] = dic['t']
                elif len(dic) == 3 and 'ang' in dic:
                    self.hop['t'][ind & (self.hop['ang'] == dic['ang'])] = dic['t']
                elif len(dic) == 3 and 'tag' in dic:
                    self.hop['t'][ind & (self.hop['tag'] == dic['tag'])] = dic['t']
                else:

                    self.hop['t'][ind & (self.hop['tag'] == dic['tag'])
                                            & (self.hop['ang'] == dic['ang'])] = dic['t']

    def get_coor_hop(self):
        '''
        Get the site coordinates in hopping space 
          only considering the nearest  neighbours hoppings.
        '''
        error_handling.empty_hop(self.hop)
        visited = np.zeros(self.lat.sites, 'u2')
        #self.lat.coor = np.sort(self.lat.coor, order=('x', 'y'))  
        self.coor_hop = np.zeros(self.lat.sites, dtype=[ ('x','f16'), ('y','f16'), ('tag', 'S1')])
        self.coor_hop['tag'] = self.lat.coor['tag']
        hop = self.hop[self.hop['n'] == 1]
        hop_down = np.copy(hop)
        hop_down['i'] = hop['j']
        hop_down['j'] = hop[ 'i']
        hop_down['ang'] = -180 + hop['ang']
        hop = np.concatenate([hop, hop_down])
        i_visit = np.min(hop['i'])
        while True:
            hs = hop[hop['i'] == i_visit]
            for h in hs:
                if visited[h['j']] == 2:
                    continue
                self.coor_hop['x'][h['j']] = self.coor_hop['x'][i_visit] + \
                    h['t'].real*cos(PI / 180 * h['ang'])
                self.coor_hop['y'][h['j']] = self.coor_hop['y'][i_visit] + \
                    h['t'].real*sin(PI / 180 * h['ang'])
                visited[h['j']] = 1
            visited[i_visit] = 2
            explored = np.argwhere(visited == 1)
            if not explored.any():
                break
            i_visit = explored[0]

    def get_ham(self):
        '''
        Get the Tight-Binding Hamiltonian.
        '''
        error_handling.empty_hop(self.hop)
        error_handling.hop_sites(self.hop, self.lat.sites)
        self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                        shape=(self.lat.sites, self.lat.sites)) \
                       + sparse.csr_matrix((self.hop['t'].conj(), (self.hop['j'], self.hop['i'])), 
                                                        shape=(self.lat.sites, self.lat.sites))
        if self.onsite.size == self.lat.sites:
            self.ham += sparse.diags(self.onsite, 0)

    def get_eig(self, eigenvec=False):
        '''
        Get the eigenergies, eigenvectors and polarisation.

        :param eigenvec: Bool. Default value False. Get the eigenvectors.
        '''
        error_handling.empty_ham(self.ham)
        error_handling.boolean(eigenvec, 'eigenvec')
        if eigenvec:
            if (self.ham.H != self.ham).nnz:
                self.en, self.vn = LA.eig(self.ham.toarray())
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
                self.vn = self.vn[:, ind]
            else:
                self.en, self.vn = LA.eigh(self.ham.toarray())
            self.intensity = np.abs(self.vn) ** 2
            self.pola = np.zeros((self.lat.sites, len(self.lat.tags)))
            for i, tag in enumerate(self.lat.tags):
                self.pola[:, i] = np.sum(self.intensity[self.lat.coor['tag'] == tag, :], axis=0)
        else:
            if (self.ham.H != self.ham).nnz:
                self.en = LA.eigvals(self.ham.toarray())
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
            else:
                self.en = LA.eigvalsh(self.ham.toarray())

    def get_ipr(self):
        '''
        Get Inverse Participation Ratio: :math:`IPR_n = |\sum_i\psi_i^{n}|^4`.

        :returns:
            * **ipr** -- IPR.
        '''
        error_handling.empty_vn(self.vn)
        self.ipr = np.sum(self.intensity ** 2, axis=0)

    def get_intensity_pola_max(self, tag_pola):
        '''
        Get the state with largest polarization on one sublattice.

        :param tag: Binary char. Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        error_handling.empty_vn(self.vn)
        error_handling.tag(tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmax(self.pola[:, i_tag])
        print('State with polarization: {:.5f}'.format(float(self.pola[ind, i_tag])))
        return self.intensity[:, ind]

    def get_intensity_pola_min(self, tag_pola):
        '''
        Get the state with smallest polarization on one sublattice.

        :param tag: Binary char. Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        error_handling.empty_vn(self.vn)
        error_handling.tag(tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmin(self.pola[:, i_tag])
        print('State with polarization: {:.5f}'.format(float(self.pola[ind, i_tag])))
        return self.intensity[:, ind]

    def get_intensity_en(self, lims):
        '''
        Get, if any, the intensity of the sum of the states 
        between *lims[0]* and *lims[1]*.

        :param lims: List, lims[0] energy min, lims[1] energy max.

        :returns:
            * **intensity** -- Sum of the intensities between (lims[0], lims[1]).
        '''
        error_handling.empty_vn(self.vn)
        error_handling.lims(lims)
        ind = np.where((self.en > lims[0]) & (self.en < lims[1]))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), lims[0], lims[1]))
        intensity = np.sum(self.intensity[:, ind], axis=1)
        return intensity
