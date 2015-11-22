import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import sys
import numpy.core.defchararray as npc

PI = np.pi


def test_hop(hop):
    if hop.size == 0:
        raise RuntimeError('\n\nRun method set_hop_uni() or set_hop() first\n')


def test_ham(ham):
    if not ham.nnz:
        raise RuntimeError('\n\nRun method get_ham() first.\n')


def test_en(en):
    if en.size == 0:
        raise RuntimeError('\n\nRun method get_eig() or first\n')


def test_vn(vn):
    if vn.size == 0:
        raise RuntimeError('\n\nRun method get_eig(eigenvec=True)  first\n')


def test_set_onsite(on):
    '''
    Check method *set_onsite*.

    :raises TypeError: Parameter *on* must be a list.
    :raises ValueError: Parameter *on* must be a container of real 
      and/or complex numbers.
    '''
    if not isinstance(on, list):
        raise TypeError('\n\nParameter *on* must be a list.\n')
    if not all([isinstance(o, (int, float, complex)) for o in on]):
        raise ValueError('\n\nParameter *on* must be a container of\
                                    real and/or complex numbers.\n')


def test_print_hop(n, n_max):
    '''
    Check method *print_vec_hop*.

    :raises TypeError: Parameter *n_max* must be an integer.
    :raises ValueError: Parameter *n_max* must be a positive integer.
    '''
    if not isinstance(n, int):
        raise TypeError('\n\nParameter *n_max* must be an integer.\n')
    if n < 1 or n > n_max-1:
        raise ValueError('\n\nParameter *n_max* must be a positive integer.\n')


def test_set_hop_uni(dict_hop, n_max):
    '''
    Check method *set_hop_uni*.

    :raises TypeError: Parameter *dict_hop* must be a dictionary.
    :raises ValueError: *key* must be a natural integer (0 < key <= nth max).
    :raises ValueError: Parameter *value* must be a number.
    '''
    if not isinstance(dict_hop, dict):
        raise TypeError('\n\nParameter *dict_ho* must be a dictionary\
                                  with key "n" and value "val".\n')
    for key, val in dict_hop.items():
        if not isinstance(key, int):
            raise ValueError('\n\n*dict_hop* keys must be integers.\n')
        if not 0 < key <= n_max:
            raise ValueError('\n\n*dict_hop* keys must be between 1 and\
                                    nth max distance between sites".\n')
        if not isinstance(val, (int, float, complex)):
            raise ValueError('\n\*dict_hop* values must be numbers".\n')

def test_set_hop(dict_hop, n_max):
    '''
    Check method *set_hop*.

    :raises TypeError: Parameter *dict_hop* must be a dictionary.
    :raises ValueError: *key* must be a natural integer (0 < key <= nth max).
    :raises ValueError: *value* must be a dictionary.
    :raises TypeError:  *dict_hop* values (*dic*) must be a dictionary.
    :raises TypeError: *key* must be a real number.
    :raises ValueError: *key* must be a positive number.
    :raises ValueError: *value* must be a number.
    '''
    if not isinstance(dict_hop, dict):
        raise TypeError('\n\nParameter *dict_hop* must be a dictionary\
                                  with key "n" and value a dictionary.\n')
    for key, dic in dict_hop.items():
        if not isinstance(key, int):
            raise ValueError('\n\nParameter *dict_hop* keys must be integers.\n')
        if not 0 < key <= n_max:
            raise ValueError('\n\n*dict_hop* key must be between 1 and\
                                    nth max distance between sites".\n')
        if not isinstance(dic, dict):
            raise TypeError('\n\n*dict_hop* values must be a dictionary.\n')                                 
        for a, t in dic.items():
            if not isinstance(a, (int, float)):
                raise ValueError('\n\nAngles must be real numbers".\n')
            if a < 0. or a > 180.:
                raise ValueError('\n\nAngles must be between 0 and 180".\n')
            if not isinstance(t, (int, float, complex)):
                raise ValueError('\n\nHoppings must be numbers".\n')

def test_set_hop_nearest(dict_hop):
    '''
    Check method *set_hop_nearest*.

    :raises TypeError: Parameter *dict_hop* must be a dictionary.
    :raises ValueError: *key* must be a natural integer (0 < key <= nth max).
    :raises ValueError: *value* must be a dictionary.
    :raises TypeError:  *dict_hop* values (*dic*) must be a dictionary.
    :raises TypeError: *key* must be a real number.
    :raises ValueError: *key* must be a positive number.
    :raises ValueError: *value* must be a number.
    '''
    if not isinstance(dict_hop, dict):
        raise TypeError('\n\nParameter *dict_hop* must be a dictionary\
                                  with key "n" and value a dictionary.\n')
    for key, val in dict_hop.items():
        if not isinstance(key, bytes):
            raise ValueError('\n\nParameter *dict_hop* keys must be hopping tags.\n')
        if len(key) != 2:
            raise ValueError('\n\nkey must be of length 2".\n')
        if not isinstance(val, (int, float)):
            raise ValueError('\n\nval must be numbers".\n')

def test_set_dimerization_def(hop, tag_hop, x_bottom_left, y_bottom_left):
    '''
    Check method *set_dimerization_def*.

    :raises TypeError: Parameter *alpha* must be a number.
    :raises ValueError: Parameter *alpha* must be a positive number.
    '''
    test_hop(hop)
    if not isinstance(tag_hop, bytes):
        raise TypeError('\n\n*tag_hop* must be a binary string of length 2.\n')
    if  len(tag_hop) != 2:
        raise ValueError('\n\n*tag_hop* must be a binary string of length 2.\n')
    if tag_hop not in hop['tag']:
        raise ValueError('\n\n*tag_hop* must be a hopping tag.\n')
    if not isinstance(x_bottom_left, (int, float)):
        raise TypeError('\n\n*x_bottom_left* must be a real number.\n')
    if not isinstance(y_bottom_left, (int, float)):
        raise TypeError('\n\n*y_bottom_left* must be a real number.\n')

def test_set_hop_disorder(hop, alpha):
    '''
    Check method *set_disorder_hop*.

    :raises TypeError: Parameter *alpha* must be a number.
    :raises ValueError: Parameter *alpha* must be a positive number.
    '''
    test_hop(hop)
    if not isinstance(alpha, (int, float, complex)):
        raise TypeError('\n\nParameter *alpha* must be a number.\n')
    if not alpha> 0:
        raise ValueError('\n\nParameter *alpha* must be positive.\n')


def test_set_on_disorder(on, alpha):
    '''
    Check method *set_disorder_hop*.

    :raises TypeError: Parameter *alpha* must be a number.
    :raises ValueError: Parameter *alpha* must be a positive number.
    '''
    if not isinstance(alpha, (int, float, complex)):
        raise TypeError('\n\nParameter *alpha* must be a number.\n')
    if not alpha> 0:
        raise ValueError('\n\nParameter *alpha* must be positive.\n')


def test_get_ham(hop, complex_transpose):
    '''
    Check method *get_ham*.

    :raises TypeError: Parameter *compl_trans* must be a bool.
    '''
    test_hop(hop)
    if not isinstance(complex_transpose, bool):
        raise TypeError('\n\nParameter *complex_transpose* must be a bool.\n')


def test_get_eig(ham, eigenvec):
    '''
    Check method *get_eig*.

    :raises TypeError: Parameter *eigenvec* must be a bool.
    '''
    test_ham(ham)
    if not isinstance(eigenvec, bool):
        raise TypeError('\n\nParameter *eigenvec* must be a bool.\n')


def test_get_state_pola(vn, tag_pola, tags):
    '''
    Check method *get_state_pola*.

    :raises TypeError: Parameter *tag_pola* must be a binary string.
    :raises ValueError: Parameter *tag_pola* is not a tag.

    '''
    test_vn(vn)
    if not isinstance(tag_pola, bytes):
        raise TypeError('\n\nParameter *tag_pola* must be a binary char.\n')
    if tag_pola not in tags:
        raise ValueError('\n\nParameter *tag_pola* is not a tag.\n')


def test_get_states_en(vn, e_min, e_max):
    '''
    Check method *get_states_en*.

    :raises TypeError: Parameter *e_min* must be a real number.
    :raises TypeError: Parameter *e_max* must be a real number.
    :raises ValueError: *e_min* must be smaller than *e_max*.
    '''
    test_vn(vn)
    if not isinstance(e_min, (int, float)):
        raise TypeError('\n\nParameter *e_min* must be a real number.\n')
    if not isinstance(e_max, (int, float)):
        raise TypeError('\n\nParameter *e_max* must be a real number.\n')
    if not e_min < e_max:
        raise ValueError('\n\n*e_min* must be smaller than *e_max*.\n')


class eigTB():
    '''
    Solve the Tight-Binding eigenvalue problem of a lattice defined 
    using the class **latticeTB**.

    :param lat: **latticeTB** class instance.
    '''
    def __init__(self, lat):
        if not lat.coor.size:
           raise RuntimeError('\n\nRun method get_lattice() of latticeTB first.\n')
        self.tags = lat.tags
        self.sites = lat.sites  # sites could be changed
        self.coor = np.copy(lat.coor)  # coor could be changed
        self.nx, self.ny = lat.nx, lat.ny
        self.vec_hop = self.get_hop()
        self.dist_uni = np.unique(self.vec_hop['d'])
        self.ang_uni = np.unique(self.vec_hop['a'])
        self.hop = np.array([], dtype=[('i', 'u4'), ('j', 'u4'), ('t', 'c16'), ('ang', 'u2'), ('tag', 'S2')])#  Hoppings
        self.onsite = np.zeros(self.sites, 'c16') #  Onsites energies
        self.ham = sparse.csr_matrix(([],([],[])), shape=(self.sites, self.sites))  # Hamiltonian
        self.en = np.array([], 'c16')  # Eigenenergies
        self.vn = np.array([], 'c16')  # Eigenvectors
        self.intensity = np.array([], 'f8')  # Intensities (|vn|**2)
        self.pola = np.array([], 'f8')  # sublattices polarisation (|vn^{(S)}|**2)
        self.alpha = 0.  # hopping disorder strength
        self.alpha_on = 0.  # onsite disorder strength
        self.params = {}

    def get_hop(self):
        '''
        Get the distances and the angles of the hoppings.
        '''
        dif_x = self.coor['x'] - self.coor['x'].reshape(self.sites, 1)
        dif_y = self.coor['y'] - self.coor['y'].reshape(self.sites, 1)
        dist = np.sqrt(dif_x ** 2 + dif_y ** 2).round(3)
        ang = (180 / PI * np.arctan2(dif_y, dif_x)).round(3)
        vec_hop = np.zeros(dist.shape, dtype=[('d', 'f8'),  ('a', 'f8')])
        vec_hop['d'] = dist
        vec_hop['a'] = ang
        return vec_hop

    def print_hop(self, n=5):
        '''
        Print the distances and the angles of all hoppings.

        :param n: Positive integer. Print the first nth hoppings
          distances and associated positive angles.
        '''
        test_print_hop(n, len(self.dist_uni)-1)
        print('Distances between sites:')
        for i, d in enumerate(self.dist_uni[1:n]):
            print(i+1, d)
            print('\twith positive angles:')
            positve_ang = self.vec_hop['a'][(self.vec_hop['d'] == d) &
                                                        (self.vec_hop['a'] >= 0.)]
            print('\t', np.unique(positve_ang))

    def set_onsite(self, on):
        '''
        Set onsite energies.

        :param on:  Array. Sublattice onsite energies.
        '''
        test_set_onsite(on)
        self.on = on
        for o, t in zip(on, self.tags):
            self.onsite[self.coor['tag'] == t] = o
           
    def set_hop_uni(self, dict_hop):
        '''
        Set uniform lattice hoppings.

        :param dict_hop: Dictionary with key {n} nth hopping
        :param dict_hop: Dictionary with key {n} nth hopping
          and value {val} of the hopping.
        '''
        test_set_hop_uni(dict_hop, len(self.dist_uni)-1)
        for key, t in dict_hop.items(): 
            ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[key]-1e-6) &
                                       (self.vec_hop['d'] < self.dist_uni[key]+ 1e-6))
            hop = np.zeros(len(ind), dtype=[('i', 'u4'), ('j', 'u4'), ('t', 'c16'), ('ang', 'u2'), ('tag', 'S2')])
            hop['i'] = ind[:, 0]
            hop['j'] = ind[:, 1]
            hop['t'] = t
            hop['tag'] = npc.add(self.coor['tag'][ind[:, 0]], self.coor['tag'][ind[:, 1]])
            self.hop = np.concatenate([self.hop, hop])

    def set_hop(self, dict_hop):
        '''
        Set non uniform lattice hoppings.

        :param dict_hop: Dictionary with key a tuple:(n, 'ang'} nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        test_set_hop(dict_hop, len(self.dist_uni)-1)
        
        for key, dic in dict_hop.items():
            ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[key]-1e-6) &
                                       (self.vec_hop['d'] < self.dist_uni[key]+1e-6))
            ind_up = ind[ind[:, 1] > ind[:, 0]]      
            dist_x = self.coor['x'][ind_up[:, 1]]-self.coor['x'][ind_up[:, 0]]
            dist_y = self.coor['y'][ind_up[:, 1]]-self.coor['y'][ind_up[:, 0]]
            ang_up = np.round(180 / PI * np.arctan2(dist_y, dist_x), 3)
            for a, t in dic.items():
                ind_tag = ind_up[ang_up == a]
                hop = np.zeros(len(ind), dtype=[('i', 'u4'), ('j', 'u4'), ('t', 'c16'), ('ang', 'u2'), ('tag', 'S2')])
                hop['i'] = ind_tag[:, 0]
                hop['j'] = ind_tag[:, 1]
                hop['ang'] = a
                hop['t'] = t
                hop['tag'] = npc.add(self.coor['tag'][ind_tag[:, 0]], self.coor['tag'][ind_tag[:, 1]])
                self.hop = np.concatenate([self.hop, hop])

    def set_hop_nearest(self, dict_hop):
        '''
        Set non uniform lattice hoppings.

        :param dict_hop: Dictionary with key a tuple:(n, 'tag'} nth hopping,
          and hopping value {val}.
        '''
        test_set_hop_nearest(dict_hop)
        ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[1]-1e-6) &
                                   (self.vec_hop['d'] < self.dist_uni[1]+1e-6))
        ind_up = ind[ind[:, 1] > ind[:, 0]]  
        self.hop = np.zeros(len(ind_up), dtype=[('i', 'u4'), ('j', 'u4'), ('t', 'c16'), ('ang', 'u2'), ('tag', 'S2')])
        self.hop['i'] = ind_up[:, 0]
        self.hop['j'] = ind_up[:, 1]
        self.hop['tag'] = npc.add(self.coor['tag'][ind_up[:, 0]], self.coor['tag'][ind_up[:, 1]])
        for s, t in dict_hop.items():
            self.hop['t'][self.hop['tag'] == s] = t

    def set_onsite_def(self, index, on_def):
        '''
        Set specific onsite energies.

        :param ind:  Array. Site indices. 
        :param on_def:  Array. Onsite energy values.
        '''
        test_set_onsite_def(index, on_def, self.sites)
        ind = np.array([index])
        on_def = np.array([on_def])
        for i, o in zip(index, on_def):
            self.onsite[i] = o

    def set_hop_def(self, index, hop_def):
        '''
        Set specific hoppings. 

        :param ind: Array. Hopping Indices, size :math:`(2 \times N_{def})`.
        :param on_def: Array. Onsite energy values.
        '''
        test_set_onsite_def(index, hop_def, self.sites)
        for i, t in zip(index, hop_def):
            cond = (self.hop['i'] == i[0]) & (self.hop['j'] == i[1])
            self.hop['t'][cond] = t
            cond = (self.hop['j'] == i[0]) & (self.hop['i'] == i[1])
            self.hop['t'][cond] = t

    def set_dimerization_def(self, tag_hop, x_bottom_left, y_bottom_left):
        '''
        Set dimerization defects.

        :param dict_hop: Dictionary with key a tuple:(n, 'ang'} nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        test_set_dimerization_def(self.hop, tag_hop, x_bottom_left, y_bottom_left)
        t1 = self.hop['t'][self.hop['tag'] == tag_hop][0]
        t2 = self.hop['t'][self.hop['tag'] == tag_hop[::-1]][0]
        pos_bool = (self.coor['x'][self.hop['i']] >= x_bottom_left) & (self.coor['y'][self.hop['i']] >= y_bottom_left)
        self.hop['t'][(pos_bool) & (self.hop['tag'] == tag_hop)] = t2
        self.hop['t'][(pos_bool) & (self.hop['tag'] == tag_hop[::-1])] = t1                             

    def  set_hop_disorder(self, alpha):
        '''
        Set a uniform hopping disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_set_hop_disorder(self.hop, alpha)
        self.hop['t'] *= 1. + alpha * rand.uniform(-0.5, 0.5, len(self.hop))
        self.alpha = alpha

    def set_on_disorder(self, alpha):
        '''
        Set a uniform onsite disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_set_on_disorder(self.onsite, alpha)
        self.onsite *= 1. + alpha * rand.uniform(-0.5, 0.5, self.sites)
        self.alpha_on = alpha

    def get_ham(self, complex_transpose=False):
        '''
        Get the Tight-Binding Hamiltonian.

        :param compl_trans: Default value False. Add complex transposed part to the Hamiltonian.
        '''
        test_get_ham(self.hop, complex_transpose)
        self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                        shape=(self.sites, self.sites)) \
                       + sparse.diags(self.onsite, 0)
        if complex_transpose:
            self.ham += sparse.csr_matrix((self.hop['t'].conj(), (self.hop['j'], self.hop['i'])), 
                                                        shape=(self.sites, self.sites))

    def get_eig(self, eigenvec=False):
        '''
        Get the eigenergies, eigenvectors and polarisations of the Tight-Binding model
        for non-Hermitian Hamiltonians.

        :param eigenvec: Default value False. Get the eigenvectors.
        '''
        test_get_eig(self.ham, eigenvec)
        if eigenvec:
            if (self.ham.H != self.ham).nnz:
                self.en, self.vn = LA.eig(self.ham.toarray())
            else:
                self.en, self.vn = LA.eigh(self.ham.toarray())
            ind = np.argsort(self.en.real)
            self.en = self.en[ind]
            self.vn = self.vn[:, ind]
            self.intensity = np.abs(self.vn) ** 2
            self.pola = np.zeros((self.sites, len(self.tags)))
            for i, t in enumerate(self.tags):
                self.pola[:, i] = np.sum(np.abs(self.vn[self.coor['tag'] == t, :]) ** 2, axis=0)
        else:
            if (self.ham.H != self.ham).nnz:
                self.en = LA.eigvals(self.ham.toarray())
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
            else:
                self.en = LA.eigvalsh(self.ham.toarray())

    def get_state_pola(self, tag_pola):
        '''
        Get the state with maximal polarization on one sublattice.

        :param tag: Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        test_get_state_pola(self.vn, tag_pola, self.tags)
        i_tag = self.tags == tag_pola
        ind = np.argmax(self.pola[:, i_tag])
        print('State with polarization:', self.pola[ind, i_tag])
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
        test_get_states_en(self.vn, e_min, e_max)
        ind = np.where((self.en > e_min) & (self.en < e_max))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), e_min, e_max))
        intensity = np.sum(np.abs(self.vn[:, ind]) ** 2, axis=1)
        return intensity
