import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import sys
import numpy.core.defchararray as npc
from math import sin, cos
from lattice import test_coor_empty, test_lat
PI = np.pi


def test_sys(sys):
    '''
    Check if other is an instance of the *system*.
    :raises TypeError: Parameter must be a instance of the class system.
    '''
    if not sys.__class__.__name__ == 'system':
        raise TypeError('\n\nParameter must be a instance of the class system.\n')


def test_coor(coor):
    '''
    Check if coordinates not empty.
    '''
    if coor.size == 0:
        raise RuntimeError('\n\nRun method get_lattice first\n')


def test_hop(hop):
    '''
    Check if hop not empty.
    '''
    if hop.size == 0:
        raise RuntimeError('\n\nRun method set_hopping first\n')


def test_ham(ham):
    '''
    Check if Hamiltonian.
    '''
    if not ham.nnz:
        raise RuntimeError('\n\nRun method get_ham first.\n')


def test_en(en):
    '''
    Check if eigenenergies.
    '''
    if en.size == 0:
        raise RuntimeError('\n\nRun method get_eig or first\n')


def test_vn(vn):
    '''
    Check if eigenvectors.
    '''
    if vn.size == 0:
        raise RuntimeError('\n\nRun method get_eig(eigenvec=True)  first\n')


def test_set_onsite(dict_onsite, tags):
    '''
    Check method *set_onsite*.

    :raises TypeError: Parameter dict_onsite must be a dictionary.
    :raises ValueError: Parameter dict_onsite keys must be a tag.
    :raises ValueError: Parameter dict_onsite values must be
      real and/or complex numbers.
    '''
    if not isinstance(dict_onsite, dict):
        raise TypeError('\n\nParameter dict_onsite must be a dictionary.\n')
    for tag, val in dict_onsite.items():
        if tag not in tags:
            raise ValueError('\n\nParameter dict_onsite keys must be a tag.\n')   
        if not isinstance(val, (int, float, complex)):
            raise ValueError('\n\nParameter dict_onsite values must be\n'\
                                       'real and/or complex numbers.\n')


def test_print_hopping(n, n_max):
    '''
    Check method *print_vec_hopping*.

    :raises TypeError: Parameter *n_max* must be an integer.
    :raises ValueError: Parameter *n_max* must be a positive integer.
      between 1 and n_max-1.
    '''
    if not isinstance(n, int):
        raise TypeError('\n\nParameter n_max must be an integer.\n')
    if n < 1 or n > n_max-1:
        raise ValueError('\n\nParameter n_max must be a positive integer'
                                    'between 1 and n_max-1.\n')

def test_set_hopping(list_hop, n_max):
    '''
    Check method *set_hopping*.

    :raises TypeError: Parameter *list_hop* must be a list.
    :raises TypeError: Parameter *list_hop* must be a list of dictionary.
    :raises KeyError: "n" and "t" must be dictionary keys.
    :raises KeyError: "tag" or "ang" must be a key.
    :raises KeyError: "tag" and "ang" must be a key.
    :raises ValueError: Dictionaries must be of length 2, 4, or 4.
    :raises ValueError: "n" must be between 1 and nmax"

    '''
    if not isinstance(list_hop, list):
        raise TypeError('\n\nParameter *list_hop* must be a list.\n')
    for dic in list_hop:
        if not isinstance(dic, dict):
            raise TypeError('\n\nParameter *list_hop* must be a list of dictionary.\n')
        if 'n' not in dic or 't' not in dic:
                raise KeyError('\n\n"n" and "t" must be dictionary keys.\n')
        if not isinstance(dic['n'], int):
            raise TypeError('\n\n"n" value must be an integer.\n')
        if not 0 < dic['n'] <= n_max:
            raise ValueError('\n\n"n" value must be between 1 and nmax".\n')
        if not isinstance(dic['t'], (int, float, complex)):
            raise TypeError('\n\n"t" value must be a real or complex number.\n')
        if len(dic) == 3:
            if not ('tag' not in dic or 'ang' not in dic):
                raise KeyError('\n\n"tag" or "ang" must be a key.\n')
        elif len(dic) == 4:
            if 'tag' not in dic and 'ang' not in dic:
                raise KeyError('\n\n"tag" or "ang" must be a key.\n')
        elif len(dic) > 4:
            raise ValueError('\n\nDictionaries must be of length 2, 3, or 4.\n')
        if 'tag' in dic:
            if not isinstance(dic['tag'], bytes):
                raise TypeError('\n\n"tag" value must be a binary string.\n')
            if len(dic['tag']) != 2:
                raise ValueError('\n\n"tag" value must be a binary string of length 2.\n')
        if 'ang' in dic:
            if not isinstance(dic['ang'], (int, float)):
                raise TypeError('\n\n"ang" value must be a real number.\n')
            if dic['ang'] < 0 or dic['ang'] >= 180:
                raise ValueError('\n\n"ang" value must be a real number'
                                          'between [0, 180].\n')

def test_rename_hop_tag(hop, list_hop):
    """
    Check method *rename_hop_tag*.

    :raises TypeError: *list_hop* must be a list.
    :raises TypeError: *list_hop* must be a list of dictionaries.
    :raises ValueError: "n" must be a key.
    :raises ValueError: "ang" must be a key.
    :raises ValueError: "tag_new" must be a key.
    """
    test_hop(hop)
    if not isinstance(list_hop, list):
        raise TypeError('\n\nnParameter *list_hop* must be a list.\n')   
    for dic in list_hop:
        if not isinstance(dic, dict):
            raise TypeError('\n\nnParameter *list_hop* must be a list of dictionaries.\n')
        if 'n' not in dic:
            raise ValueError('\n\n"n" must be a key.\n')
        if 'ang' not in dic:
            raise ValueError('\n\n"ang" must be a key.\n')
        if 'tag_new' not in dic:
            raise ValueError('\n\n"tag_new" must be a key.\n')


def test_change_hopping(hop, dict_hop, x_bottom_left, y_bottom_left):
    '''
    Check method *set_defect_dimer*.

    :raises TypeError: Parameter *alpha* must be a number.
    '''
    test_hop(hop)
    if not isinstance(dict_hop, dict):
        raise TypeError('\n\nParameter *dict_hop* must be a dictionary\
                                  with binary string keys and hopping values.\n')
    for key, val in dict_hop.items():
        if not isinstance(key, bytes):
            raise TypeError('\n\nDictionary key must be a binary string of length 2.\n')
        if  len(key) != 2:
            raise ValueError('\n\nDictionary key must be a binary string of length 2.\n')
        if key not in hop['tag']:
            raise ValueError('\n\nDictionary key must be a hopping tag.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\nDictionary value must be a number.\n')            
    if not isinstance(x_bottom_left, (int, float)):
        raise TypeError('\n\n*x_bottom_left* must be a real number.\n')
    if not isinstance(y_bottom_left, (int, float)):
        raise TypeError('\n\n*y_bottom_left* must be a real number.\n')


def test_set_hopping_def(hop, hopping_def, sites):
    '''
    Check method *test_set_hop_def*.

    :raises TypeError: Parameter *hopping_def* must be a dictionary
    :raises TypeError: *hopping_def* keys must be lists.
    :raises ValueError: *hopping_def* keys must be lists of length 2.
    :raises ValueError: *hopping_def* keys must be lists of integers.
    :raises TypeError: *hopping_def* keys must be lists.
    :raises ValueError: *hopping_def* keys must be integers between 0 and sites-1.
    :raises ValueError: *hopping_def* keys must be different integers between 0 and sites-1.
    :raises TypeError: *hopping_def* values must be numbers.
    '''
    if not isinstance(hopping_def, dict):
        raise TypeError('\n\nParameter hopping_def must be a dictionary.\n')
    for key, val in hopping_def.items():
        if not isinstance(key, tuple):
            raise TypeError('\n\nhopping_def keys must be lists.\n')
        if len(key) != 2:
            raise TypeError('\n\nhopping_def keys must be lists of length 2.\n')
        if not isinstance(key[0], int) or not isinstance(key[1], int):
            raise ValueError('\n\nhopping_def keys must be lists of integers.\n')
        if key[0] < 0 or key[1] < 0 or key[0] > sites-1 or key[1] > sites-1:
            raise ValueError('\n\nhopping_def keys must be integers between 0 and sites-1.\n')
        if key[0] == key[1]:
            raise ValueError('\n\nhopping_def keys must be different integers between 0 and sites-1.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\nhopping_def values must be numbers.\n')


def test_set_onsite_def(onsite_def, sites):
    '''
    Check method *test_set_ons_def*.

    :raises TypeError: Parameter *onsite_def* must be a dictionary.
    :raises TypeError: *onsite_def* keys must be integers.
    :raises TypeError: *onsite_def* values must be numbers.
    :raises ValueError: *onsite_def* keys must be integers between :math:`[0, sites)`.
    '''
    if not isinstance(onsite_def, dict):
        raise TypeError('\n\nParameter onsite_def must be a dictionary.\n')
    for key, val in onsite_def.items():
        if not isinstance(key, int):
            raise TypeError('\n\nonsite_def keys must be integers.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\nonsite_def values must be numbers.\n')
        if key < 0 or key > sites-1:
            raise ValueError('\n\nonsite_def keys must be integers between 0 and sites-1.\n')


def test_alpha(alpha):
    '''
    Check parameter *alpha*.

    :raises TypeError: Parameter *alpha* must be a number.
    '''
    if not isinstance(alpha, (int, float, complex)):
        raise TypeError('\n\nParameter *alpha* must be a number.\n')


def test_get_coor_hop(hop):
    '''
    Check method *get_ham*.

    :raises RunTimeError: 'Works only for nearest hoppings'.
    '''
    if len(hop['n'] == 1) == 0:
        raise ValueError('\n\nWorks only for nearest hoppings.\n')


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
    if not isinstance(tag_pola, bytes):
        raise TypeError('\n\nParameter *tag_pola* must be a binary char.\n')
    if tag_pola not in tags:
        raise ValueError('\n\nParameter *tag_pola* is not a tag.\n')


def test_en_lims(vn, en_lims):
    '''
    Check method *get_intensity_en*.

    :raises TypeError: Parameter en_lims must be a list.
    :raises TypeError: Parameter *en_lims[0]* must be a real number.
    :raises TypeError: Parameter *en_lims[1]* must be a real number.
    :raises ValueError: *en_lims* must be a list of length 2.
    :raises ValueError: *e_min* must be smaller than *e_max*.
    '''
    if not isinstance(en_lims, list):
        raise TypeError('\n\nParameter en_lims must be a list.\n')
    if not isinstance(en_lims[0], (int, float)):
        raise TypeError('\n\nParameter *en_lims[0]* must be a real number.\n')
    if not isinstance(en_lims[1], (int, float)):
        raise TypeError('\n\nParameter *en_lims[1]* must be a real number.\n')
    if len(en_lims) != 2:
        raise ValueError('\n\n*en_lims* must be a list of length 2.\n')
    if not en_lims[0] < en_lims[1]:
        raise ValueError('\n\n*en_lims[0]* must be smaller than *en_lims[1]*.\n')


class system():
    '''
    Solve the Tight-Binding eigenvalue problem of a lattice defined 
    using the class **lattice**.

    :param lat: **lattice** class instance.
    '''
    def __init__(self, lat):
        test_lat(lat)
        test_coor_empty(lat.coor)
        self.lat = lat
        self.vec_hop = self.get_hop()
        self.dist_uni = np.unique(self.vec_hop['d'])
        self.ang_uni = np.unique(self.vec_hop['a'])
        self.coor_hop = np.array([], dtype=[ ('x','f16'), ('y','f16'), ('tag','S1')])
        self.hop = np.array([], dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                    ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')]) #  Hoppings
        self.onsite = np.zeros(self.lat.sites, 'c16') #  Onsite energies
        self.ham = sparse.csr_matrix(([],([],[])), shape=(self.lat.sites, self.lat.sites))  # Hamiltonian
        self.en = np.array([], 'c16')  # Eigenenergies
        self.vn = np.array([], 'c16')  # Eigenvectors
        self.intensity = np.array([], 'f8')  # Intensities (|vn|**2)
        self.pola = np.array([], 'f8')  # sublattices polarisation (|vn^{(S)}|**2)
        self.alpha = 0.  # hopping disorder strength
        self.alpha_onsite = 0.  # onsite disorder strength
        self.params = {}

    def get_hop(self):
        '''
        Get the distances and the angles of the hoppings.
        '''
        dif_x = self.lat.coor['x'] - self.lat.coor['x'].reshape(self.lat.sites, 1)
        dif_y = self.lat.coor['y'] - self.lat.coor['y'].reshape(self.lat.sites, 1)
        dist = np.sqrt(dif_x ** 2 + dif_y ** 2).round(3)
        ang = (180 / PI * np.arctan2(dif_y, dif_x)).round(3)
        vec_hop = np.zeros(dist.shape, dtype=[('d', 'f16'),  ('a', 'f16')])
        vec_hop['d'] = dist
        vec_hop['a'] = ang
        return vec_hop

    def print_hopping(self, n=5):
        '''
        Print the distances and the angles of all hoppings.

        :param n: Positive integer. Print the first nth hopping
          distances and associated positive angles.
        '''
        n_max = len(self.dist_uni) - 1
        test_print_hopping(n, n_max)
        print('\n{} different distances between nodes:'.format(n_max))
        print('\nDistances between sites:')
        if n > n_max:
            n = n_max
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
            positve_ang = self.vec_hop['a'][(self.vec_hop['d'] == d) &
                                                              (self.vec_hop['a'] >= 0.) &
                                                              (self.vec_hop['a'] < 180.)]
            print('\t', np.unique(positve_ang))

    def set_onsite(self, dict_onsite):
        '''
        Set onsite energies.

        :param on:  Array. Sublattice onsite energies.
        '''
        test_set_onsite(dict_onsite, self.lat.tags)
        for tag, on in dict_onsite.items():
            self.onsite[self.lat.coor['tag'] ==tag] = on

    def set_hopping(self, list_hop):
        '''
        Set lattice hoppings.

        :param list_hop: List of dictionaries, Dictionary with key a tuple:(n, 'ang') nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        test_set_hopping(list_hop, len(self.dist_uni) - 1)
        list_n = np.unique([dic['n'] for dic in list_hop])
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

    def rename_hop_tag(self, list_hop):
        test_rename_hop_tag(self.hop, list_hop)
        print(list_hop)
        for dic in list_hop:
            self.hop['tag'][(self.hop['n'] == dic['n']) & (self.hop['tag'] == dic['tag']) & 
                                  (self.hop['ang'] == dic['ang'])] = dic['tag_new']

    def set_hopping_tag(self, list_hop):
        for dic in list_hop:
            self.hop['t'][self.hop['tag'] == dic['tag']] = dic['t']

    def set_onsite_def(self, onsite_def):
        '''
        Set specific onsite energies.

        :param dict_ons_def:  Dictionary. key: site indices, val: onsite values. 
        '''
        test_set_onsite_def(onsite_def, self.lat.sites)
        for i, o in onsite_def.items():
            self.onsite[i] = o

    def set_hopping_def(self, hopping_def):
        '''
        Set specific hoppings. 

        :param dict_hop_def:  Dictionary. key: hopping indices, val: hopping values. 
        '''
        test_hop(self.hop)
        test_set_hopping_def(self.hop, hopping_def, self.lat.sites)

        for key, val in hopping_def.items():
            cond = (self.hop['i'] == key[0]) & (self.hop['j'] == key[1])
            self.hop['t'][cond] = val
            cond = (self.hop['j'] == key[0]) & (self.hop['i'] == key[1])
            self.hop['t'][cond] = val

    def  set_hopping_dis(self, alpha):
        '''
        Set a uniform hopping disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_hop(self.hop)
        test_alpha(alpha)
        self.hop['t'] *= 1. + alpha * rand.uniform(-1., 1., len(self.hop))
        self.alpha = alpha

    def set_onsite_dis(self, alpha):
        '''
        Set a uniform onsite disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_alpha(alpha)
        self.onsite *= 1. + alpha * rand.uniform(-1., 1., self.lat.sites)
        self.alpha_onsite = alpha

    def set_hopping_new(self, list_hop, x_bottom_left=0, y_bottom_left=0):
        '''
        Set a dimerization defect.

        :param dict_hop: Dictionary with key a tuple:(n, 'ang'} nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        #test_set_defect_dim(self.hop, dict_hop, x_bottom_left, y_bottom_left)
        for dic in list_hop:
            ind = (self.lat.coor['x'][self.hop['i']] >= x_bottom_left) & \
                    (self.lat.coor['y'][self.hop['i']] >= y_bottom_left)
            self.hop['t'][ind & (self.hop['tag'] == dic['tag'])] = dic['t']

    def get_coor_hop(self):
        '''
        Get the site coordinates in hopping space 
          only considering the nearest  neighbours hoppings.
        '''
        test_hop(self.hop)
        test_get_coor_hop(self.hop)
        visited = np.zeros(self.lat.sites, 'u2')
        self.lat.coor = np.sort(self.lat.coor, order=('x', 'y'))  
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
        test_hop(self.hop)
        self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                        shape=(self.lat.sites, self.lat.sites)) \
                       + sparse.csr_matrix((self.hop['t'].conj(), (self.hop['j'], self.hop['i'])), 
                                                        shape=(self.lat.sites, self.lat.sites)) \
                       + sparse.diags(self.onsite, 0)

    def get_eig(self, eigenvec=False):
        '''
        Get the eigenergies, eigenvectors and polarisationsite of the Tight-Binding model
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
            tags_uni = np.unique([dic['tag'] for dic in self.lat.unit_cell])
            self.pola = np.zeros((self.lat.sites, len(tags_uni)))
            for i, tag in enumerate(tags_uni):
                self.pola[:, i] = np.sum(np.abs(self.vn[self.lat.coor['tag'] == tag, :]) ** 2, axis=0)
        else:
            if (self.ham.H != self.ham).nnz:
                self.en = LA.eigvals(self.ham.toarray())
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
            else:
                self.en = LA.eigvalsh(self.ham.toarray())

    def get_intensity_pola(self, tag_pola):
        '''
        Get the state with maximal polarization on one sublattice.

        :param tag: Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        test_vn(vn)
        test_get_state_pola(tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmax(self.pola[:, i_tag])
        print('State with polarization:', self.pola[ind, i_tag])
        return self.intensity[:, ind]

    def get_intensity_en(self, en_lims):
        '''
        Get, if any, the intensity of the sum of the states 
        between *en_lims[0]* and *en_lims[1]*.

        :param en_lims: List, en_lims[0] energy min, en_lims[1] energy max.

        :returns:
            * **intensity** -- Sum of the intensities between *e_min* and *e_max*.
        '''
        test_vn(vn)
        test_en_lims(en_lims)
        ind = np.where((self.en > en_lims[0]) & (self.en < en_lims[1]))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), en_lims[0], en_lims[1]))
        intensity = np.sum(np.abs(self.vn[:, ind]) ** 2, axis=1)
        return intensity





        '''
        for dic in list_hop:
            ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[dic_n['n']] - 1e-4) &
                                       (self.vec_hop['d'] < self.dist_uni[dic_n['n']] + 1e-4))
            ind_up = ind[ind[:, 1] > ind[:, 0]]
            hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                               ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
            hop['n'] = dic_n['n']
            hop['i'] = ind_up[:, 0]
            hop['j'] = ind_up[:, 1]
            hop['ang'] = self.vec_hop['a'][ind_up[:, 0], ind_up[:, 1]]
            hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], self.lat.coor['tag'][ind_up[:, 1]])
            for dic_hop in  dic['hop']:
                hop['t'][hop['tag'] == dic_hop['tag']] = dic_hop['t']
            self.hop = np.concatenate([self.hop, hop])
        '''