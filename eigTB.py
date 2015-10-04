import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import sys


class eigTB():
    '''
    Solve the Tight-Binding eigenvalue problem of a lattice defined 
    using the class latticeTB.

    :param coor: Sites coordinates (class **latticeTB** output).
    :param tags: Sites tags (class **latticeTB** output).
    '''
    def __init__(self, lat):
        if not (lat.coor['x'].any() or lat.coor['y'].any()):
           raise Exception('\n\nRun method get_lattice() from latticeTB first.\n')
        self.tags = lat.tags
        self.sites = np.copy(lat.sites)  # sites could be changed
        self.coor = np.copy(lat.coor)  # coor could be changed
        self.on = np.array([], 'c16') #  Sublattice onsite energies
        self.hop = np.array([], {'names': ['i', 'j', 't'], 'formats': ['u4', 'u4', 'c16']}) #  Hoppings
        self.onsite = np.zeros(self.sites, 'c16') #  Onsites energies
        self.ham = sparse.csr_matrix(([],([],[])), shape=(self.sites, self.sites))  # Hamiltonian
        self.en = np.array([], 'c16')  # Eigenenergies
        self.vn = np.array([], 'c16')  # Eigenvectors
        self.intensity = np.array([], 'f8')  # Intensities (|vn|**2)
        self.pola = np.array([], 'f8')  # sublattices polarisation (|vn^{(S)}|**2)
        self.alpha_hop = 0.  # hopping disorder strength
        self.alpha_on = 0.  # onsite disorder strength
        self.params = {}

    def set_onsite(self, on):
        '''
        Set onsite energies.

        :param on:  Array. Sublattices onsite energies.
        '''
        self.on = on
        for o, t in zip(on, self.tags):
            self.onsite[self.coor['tag'] == t] = o

    def set_onsite_def(self, ind, on_def):
        '''
        Set specific onsite energies.

        :param ind:  Array. Site indices. 
        :param on_def:  Array. Onsite energy values.
        '''
        ind = np.array([ind])
        on_def = np.array([on_def])
        for i, o in zip(ind, on_def):
            self.onsite[i] = o

    def set_hop_def(self, ind, hop_def):
        r'''
        Set specific hoppings. 

        :param ind: Array. Hopping Indices, size :math:`(2 \times N_{def})`.
        :param on_def: Array. Onsite energy values.
        '''
        if not self.hop['t'].any():
            raise Exception('\n\nRun method set_hop() first.\n')
        for i, t in zip(ind, hop_def):
            cond = (self.hop['i'] == i[0]) & (self.hop['j'] == i[1])
            self.hop['t'][cond] = t

    def set_disorder_hop(self, alpha=0):
        '''
        Set a generic disorder. 

        :param alpha: Default value 0. Stength of the disorder.
        '''
        if not self.hop['t'].any():
            raise Exception('\n\nRun method set_hop() first.\n')
        self.hop['t'] *= 1+ alpha * rand.uniform(-0.5, 0.5, len(self.hop['tag']))
        self.alpha_hop = alpha
        print('dis',  self.alpha_hop)

    def set_disorder_on(self, alpha=0):
        '''
        Set a generic disorder. 

        :param alpha: Default value 0. Stength of the disorder.
        '''
        self.onsite *= 1+ alpha * rand.uniform(-0.5, 0.5, len(self.sites))
        self.alpha_on = alpha

    def get_ham(self):
        '''
        Get the Tight-Binding Hamiltonian.
        '''
        if not self.hop['t'].any():
            raise Exception('\n\nRun method set_hop() first.\n')
        self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                        shape=(self.sites, self.sites)) \
                       + sparse.csr_matrix((self.hop['t'], (self.hop['j'], self.hop['i'])), 
                                                        shape=(self.sites, self.sites)).conj() \
                       + sparse.diags(self.onsite, 0)

    def get_eig(self):
        '''
        Get the eigenergies, eigenvectors and polarisations of the Tight-Binding model
        for non-Hermitian Hamiltonians.
        '''
        if not self.ham.nnz:
            raise Exception('\n\nRun method get_ham() first.\n')
        self.en, self.vn = LA.eig(self.ham.toarray())
        ind = np.argsort(self.en.real)
        self.en = self.en[ind]
        self.vn = self.vn[:, ind]
        self.intensity = np.abs(self.vn) ** 2
        self.pola = np.zeros((self.sites, len(self.tags)))
        for i, t in enumerate(self.tags):
            self.pola[:, i] = np.sum(np.abs(self.vn[self.coor['tag'] == t, :]) ** 2, axis=0)

    def get_eigh(self):
        '''
        Get the eigenergies, eigenvectors and polarisations of the Tight-Binding model
        for Hermitian Hamiltonians.
        '''
        if not self.ham.nnz:
            raise Exception('\n\nRun method get_ham() first.\n')
        self.en, self.vn = LA.eigh(self.ham.toarray())
        self.intensity = np.abs(self.vn) ** 2
        self.pola = np.zeros((self.sites, len(self.tags)))
        for i, t in enumerate(self.tags):
            self.pola[:, i] = np.sum(np.abs(self.vn[self.coor['tag'] == t, :]) ** 2, axis=0)

    def get_state_pola(self, tag):
        '''
        Get the state with maximal polarization on one sublattice.

        :param tag: Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        if not self.pola.any():
            raise Exception('\n\nRun method get_eig() or get_eigh() first.\n')
        ind = np.argmax(self.pola[:, self.tags == tag])
        print('State with polarization:', self.pola[ind, self.tags == tag])
        return self.intensity[:, ind]

    def get_states_en(self, e_min, e_max):
        '''
        Get, if any, the intensity of the sum of the states 
        between *e_min* and *e_max*.

        :param e_min: Energy min.
        :param e_max: Energy max.

        :returns:
            * **intensity** -- Sum of the intensities between *e_min* and *e_max*.
        '''
        if not self.en.any():
            raise Exception('\n\nRun method get_eig() or get_eigh() first.\n')
        ind = np.where((self.en > e_min) & (self.en < e_max))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), e_min, e_max))
        intensity = np.sum(np.abs(self.vn[:, ind]) ** 2, axis=1)
        return intensity
