import numpy as np
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import sys
import numpy.core.defchararray as npc
from math import sin, cos
PI = np.pi


def test_hop(hop):
    '''
    Check if hoppings.
    '''
    if hop.size == 0:
        raise RuntimeError('\n\nRun method set_hop_uni() or set_hop() first\n')


def test_ham(ham):
    '''
    Check if Hamiltonian.
    '''
    if not ham.nnz:
        raise RuntimeError('\n\nRun method get_ham() first.\n')


def test_en(en):
    '''
    Check if eigenenergies.
    '''
    if en.size == 0:
        raise RuntimeError('\n\nRun method get_eig() or first\n')


def test_vn(vn):
    '''
    Check if eigenvectors.
    '''
    if vn.size == 0:
        raise RuntimeError('\n\nRun method get_eig(eigenvec=True)  first\n')


def test_set_ons(on):
    '''
    Check method *set_ons*.

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
        raise TypeError('\n\nParameter *dict_hop* must be a dictionary\
                                  with key "n" and value "val".\n')
    for key, val in dict_hop.items():
        if not isinstance(key, int):
            raise ValueError('\n\n*dict_hop* keys must be integers.\n')
        if not 0 < key <= n_max:
            raise ValueError('\n\n*dict_hop* keys must be between 1 and nth max".\n')
        if not isinstance(val, (int, float, complex)):
            raise ValueError('\n\*dict_hop* values must be numbers".\n')

def test_set_hop(list_hop, n_max):
    '''
    Check method *set_hop*.

    :raises TypeError: Parameter *list_hop* must be a list.
    :raises TypeError: Parameter *list_hop* must be a list of dictionary.
    :raises TypeError: "n" must be a key.
    :raises TypeError: "hop" must be a key.
    :raises ValueError: "tag" must be a key.
    :raises ValueError: "t" must be a key.

    '''
    if not isinstance(list_hop, list):
        raise TypeError('\n\nParameter *list_hop* must be a list.\n')
    for dic in list_hop:
        if not isinstance(dic, dict):
            raise TypeError('\n\nParameter *list_hop* must be a list of dictionary.\n')
        if 'n' not in dic:
                raise ValueError('\n\n"n" must be a key.\n')
        if 'hop' not in dic:
                raise ValueError('\n\n"hop" must be a key.\n')
        if not 0 < dic['n'] <= n_max:
            raise ValueError('\n\n"n" must be between 1 and\
                                    nth max".\n')
        for d in dic['hop']:
            if 't' not in d:
                raise ValueError('\n\n"n" must be a key.\n')
            if 'tag' not in d:
                raise ValueError('\n\n"t" must be a key.\n')


def test_rename_hop_tag(hop, list_hop):
    """
    Check method *rename_hop_tag*.

    :raises TypeError: *list_hop* must be a list.
    :raises TypeError: *list_hop* must be a list of dictionaries.
    :raises ValueError: "n" must be a key.
    :raises ValueError: "ang" must be a key.
    :raises ValueError: "tag" must be a key.
    :raises ValueError: "new_tag" must be a key.
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
        if 'tag' not in dic:
            raise ValueError('\n\n"tag" must be a key.\n')
        if 'new_tag' not in dic:
            raise ValueError('\n\n"new_tag" must be a key.\n')


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

def test_set_dimerization_def(hop, dict_hop, x_bottom_left, y_bottom_left):
    '''
    Check method *set_dimerization_def*.

    :raises TypeError: Parameter *alpha* must be a number.
    :raises ValueError: Parameter *alpha* must be a positive number.
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


def test_set_hop_def(hop, dict_hop, sites):
    '''
    Check method *test_set_hop_def*.

    :raises TypeError: Parameter *dict_hop* must be a dictionary
    :raises TypeError: *dict_hop* keys must be lists.
    :raises ValueError: *dict_hop* keys must be lists of length 2.
    :raises ValueError: *dict_hop* keys must be lists of integers.
    :raises TypeError: *dict_hop* keys must be lists.
    :raises ValueError: *dict_hop* keys must be integers between 0 and sites-1.
    :raises ValueError: *dict_hop* keys must be different integers between 0 and sites-1.
    :raises TypeError: *dict_hop* values must be numbers.
    '''
    test_hop(hop)
    if not isinstance(dict_hop, dict):
        raise TypeError('\n\nParameter *dict_hop* must be a dictionary.\n')
    for key, val in dict_hop.items():
        if not isinstance(key, tuple):
            raise TypeError('\n\n*dict_hop* keys must be lists.\n')
        if len(key) != 2:
            raise TypeError('\n\n*dict_hop* keys must be lists of length 2.\n')
        if not isinstance(key[0], int) or not isinstance(key[1], int):
            raise ValueError('\n\n*dict_hop* keys must be lists of integers.\n')
        if key[0] < 0 or key[1] < 0 or key[0] > sites-1 or key[1] > sites-1:
            raise ValueError('\n\n*dict_hop* keys must be integers between 0 and sites-1.\n')
        if key[0] == key[1]:
            raise ValueError('\n\n*dict_hop* keys must be different integers between 0 and sites-1.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\n*dict_hop* values must be numbers.\n')


def test_set_ons_def(dict_ons, sites):
    '''
    Check method *test_set_ons_def*.

    :raises TypeError: Parameter *dict_ons* must be a dictionary.
    :raises TypeError: *dict_ons* keys must be integers.
    :raises TypeError: *dict_ons* keys must be numbers.
    :raises ValueError: *dict_ons* keys must be integers between 0 and sites-1.
    '''
    if not isinstance(dict_ons, dict):
        raise TypeError('\n\nParameter *dict_ons* must be a dictionary.\n')
    for key, val in dict_ons.items():
        if not isinstance(key, int):
            raise TypeError('\n\n*dict_ons* keys must be integers.\n')
        if not isinstance(val, (int, float, complex)):
            raise TypeError('\n\n*dict_ons* values must be numbers.\n')
        if key < 0 or key > sites-1:
            raise ValueError('\n\n*dict_ons* keys must be integers between 0 and sites-1.\n')


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


def test_set_ons_disorder(on, alpha):
    '''
    Check method *set_disorder_hop*.

    :raises TypeError: Parameter *alpha* must be a number.
    :raises ValueError: Parameter *alpha* must be a positive number.
    '''
    if not isinstance(alpha, (int, float, complex)):
        raise TypeError('\n\nParameter *alpha* must be a number.\n')
    if not alpha> 0:
        raise ValueError('\n\nParameter *alpha* must be positive.\n')


def test_get_coor_hop(hop):
    '''
    Check method *get_ham*.

    :raises RunTimeError: 'Works only for nearest hoppings'.
    '''
    test_hop(hop)
    if len(hop['n'] == 1) == 0:
        raise RuntimeError('\n\nWorks only for nearest hoppings.\n')

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


def test_get_states_en(vn, en_lims):
    '''
    Check method *get_states_en*.

    :raises TypeError: Parameter en_lims must be a list.
    :raises TypeError: Parameter *en_lims[0]* must be a real number.
    :raises TypeError: Parameter *en_lims[1]* must be a real number.
    :raises ValueError: *e_min* must be smaller than *e_max*.
    '''
    test_vn(vn)
    if not isinstance(en_lims, list):
        raise TypeError('\n\nParameter en_lims must be a list.\n')
    if not isinstance(en_lims[0], (int, float)):
        raise TypeError('\n\nParameter *en_lims[0]* must be a real number.\n')
    if not isinstance(en_lims[1], (int, float)):
        raise TypeError('\n\nParameter *en_lims[1]* must be a real number.\n')
    if not en_lims[0] < en_lims[1]:
        raise ValueError('\n\n*en_lims[0]* must be smaller than *en_lims[1]*.\n')


class eigTB():
    '''
    Solve the Tight-Binding eigenvalue problem of a lattice defined 
    using the class **latticeTB**.

    :param lat: **latticeTB** class instance.
    '''
    def __init__(self, lat):
        if not lat.coor.size:
           raise RuntimeError('\n\nRun method get_lattice() of latticeTB first.\n')
        self.lat = lat
        self.vec_hop = self.get_hop()
        self.dist_uni = np.unique(self.vec_hop['d'])
        self.ang_uni = np.unique(self.vec_hop['a'])
        self.coor_hop = np.array([], dtype=[ ('x','f16'), ('y','f16'), ('tag','S1')])
        self.hop = np.array([], dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                    ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')]) #  Hoppings
        self.ons = np.zeros(self.lat.sites, 'c16') #  Onsite energies
        self.ham = sparse.csr_matrix(([],([],[])), shape=(self.lat.sites, self.lat.sites))  # Hamiltonian
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
        dif_x = self.lat.coor['x'] - self.lat.coor['x'].reshape(self.lat.sites, 1)
        dif_y = self.lat.coor['y'] - self.lat.coor['y'].reshape(self.lat.sites, 1)
        dist = np.sqrt(dif_x ** 2 + dif_y ** 2).round(3)
        ang = (180 / PI * np.arctan2(dif_y, dif_x)).round(3)
        vec_hop = np.zeros(dist.shape, dtype=[('d', 'f16'),  ('a', 'f16')])
        vec_hop['d'] = dist
        vec_hop['a'] = ang
        return vec_hop

    def print_hop(self, n=5):
        '''
        Print the distances and the angles of all hoppings.

        :param n: Positive integer. Print the first nth hopping
          distances and associated positive angles.
        '''
        test_print_hop(n, len(self.dist_uni)-1)
        print('Distances between sites:')
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
                                                        (self.vec_hop['a'] >= 0.)]
            print('\t', np.unique(positve_ang))

    def set_ons(self, on):
        '''
        Set onsite energies.

        :param on:  Array. Sublattice onsite energies.
        '''
        test_set_ons(on)
        self.on = on
        for o, t in zip(on, self.lat.tags):
            self.ons[self.lat.coor['tag'] == t] = o

    def set_hop_uni(self, list_hop):
        '''
        Set uniform lattice hoppings.

        :param list_hop: Lisi of dictionaries.
           Dictionary with keys 'n' and 't' nth hopping and hopping value.
        '''
        #test_set_hop_uni(dict_hop, len(self.dist_uni)-1)
        for dic in list_hop: 
            ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[dic['n']]-1e-4) &
                                       (self.vec_hop['d'] < self.dist_uni[dic['n']]+ 1e-4))
            hop = np.zeros(len(ind), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                           ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
            hop['n'] = dic['n']
            hop['i'] = ind[:, 0]
            hop['j'] = ind[:, 1]
            hop['t'] = dic['t']
            hop['ang'] = self.vec_hop['a'][ind[:, 0], ind[:, 1]]
            hop['tag'] = npc.add(self.lat.coor['tag'][ind[:, 0]], self.lat.coor['tag'][ind[:, 1]])
            self.hop = np.concatenate([self.hop, hop])

    def set_hop(self, list_hop):
        '''
        Set non uniform lattice hoppings.

        :param list_hop: List of dictionaries, Dictionary with key a tuple:(n, 'ang') nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        test_set_hop(list_hop, len(self.dist_uni) - 1)
        for dic_n in list_hop:
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
            for dic_hop in  dic_n['hop']:
                hop['t'][hop['tag'] == dic_hop['tag']] = dic_hop['t']
            self.hop = np.concatenate([self.hop, hop])

    def rename_hop_tag(self, list_hop):
        for dic in list_hop:
            self.hop['tag'][(self.hop['n'] == dic['n']) & (self.hop['ang'] == dic['ang'])] = dic['tag_new']
            self.hop['tag'][(self.hop['n'] == dic['n']) & (self.hop['ang'] == -180 + dic['ang'])] = dic['tag_new']

    def set_hop_with_tag(self, list_hop):
        for dic in list_hop:
            self.hop['t'][self.hop['tag'] == dic['tag']] = dic['t']

    def set_hop_nearest(self, dict_hop):
        '''
        Set only nearest hoppings via hopping tags.

        :param dict_hop: Dictionary with key a tuple:(n, 'tag'} nth hopping,
          and hopping value {val}.
        '''
        test_set_hop_nearest(dict_hop)
        ind = np.argwhere((self.vec_hop['d'] > self.dist_uni[1]-1e-4) &
                                   (self.vec_hop['d'] < self.dist_uni[1]+1e-4))
        ind_up = ind[ind[:, 1] > ind[:, 0]]  
        self.hop = np.zeros(len(ind_up), dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                                                  ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
        self.hop['n'] = 1
        self.hop['i'] = ind_up[:, 0]
        self.hop['j'] = ind_up[:, 1]
        self.hop['ang'] = self.vec_hop['a'][ind_up[:, 0], ind_up[:, 1]]
        self.hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], self.lat.coor['tag'][ind_up[:, 1]])
        for s, t in dict_hop.items():
            self.hop['t'][self.hop['tag'] == s] = t

    def set_ons_def(self, dict_ons_def):
        '''
        Set specific onsite energies.

        :param dict_ons_def:  Dictionary. key: site indices, val: onsite values. 
        '''
        test_set_ons_def(dict_ons_def, self.lat.sites)
        for i, o in dict_ons_def.items():
            self.ons[i] = o

    def set_hop_def(self,dict_hop_def):
        '''
        Set specific hoppings. 

        :param dict_hop_def:  Dictionary. key: hopping indices, val: hopping values. 
        '''
        test_set_hop_def(self.hop, dict_hop_def, self.lat.sites)
        for key, val in dict_hop_def.items():
            cond = (self.hop['i'] == key[0]) & (self.hop['j'] == key[1])
            self.hop['t'][cond] = val
            cond = (self.hop['j'] == key[0]) & (self.hop['i'] == key[1])
            self.hop['t'][cond] = val

    def set_dimerization_def(self, dict_hop, x_bottom_left, y_bottom_left):
        '''
        Set a dimerization defect.

        :param dict_hop: Dictionary with key a tuple:(n, 'ang'} nth hopping,
          associated positive angle, and hopping value {val}.
        '''
        test_set_dimerization_def(self.hop, dict_hop, x_bottom_left, y_bottom_left)
        for key, val in dict_hop.items():
            ind = (self.lat.coor['x'][self.hop['i']] >= x_bottom_left) & \
                    (self.lat.coor['y'][self.hop['i']] >= y_bottom_left)
            self.hop['t'][ind & (self.hop['tag'] == key)] = val

    def  set_hop_disorder(self, alpha):
        '''
        Set a uniform hopping disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_set_hop_disorder(self.hop, alpha)
        self.hop['t'] *= 1. + alpha * rand.uniform(-0.5, 0.5, len(self.hop))
        self.alpha = alpha

    def set_ons_disorder(self, alpha):
        '''
        Set a uniform onsite disorder. 

        :param alpha: Stength of the disorder.
        '''
        test_set_ons_disorder(self.ons, alpha)
        self.ons *= 1. + alpha * rand.uniform(-0.5, 0.5, self.lat.sites)
        self.alpha_on = alpha

    def get_coor_hop(self):
        '''
        Get the site coordinates in hopping space 
          considering just the nearest  hoppings.

        :param n: Unsigned int. Specific hopping.
        '''
        test_get_coor_hop(self.hop)
        visited = np.zeros(self.lat.sites, 'u2')
        self.coor_hop = np.zeros(self.lat.sites, dtype=[ ('x','f16'), ('y','f16'), ('tag', 'S1')])
        self.coor_hop['tag'] = self.lat.coor['tag']
        hop = self.hop[self.hop['n'] == 1]
        hop_down = np.copy(hop)
        hop_down['i'] = hop['j']
        hop_down['j'] = hop['i']
        hop_down['ang'] = -180 + hop['ang']
        hop = np.concatenate([hop, hop_down])
        #print(hop)
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

    def get_ham(self, complex_transpose=False):
        '''
        Get the Tight-Binding Hamiltonian.

        :param compl_trans: Default value False. Add complex transposed part to the Hamiltonian.
        '''
        test_get_ham(self.hop, complex_transpose)
        self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                        shape=(self.lat.sites, self.lat.sites)) \
                       + sparse.diags(self.ons, 0)
        if complex_transpose:
            self.ham += sparse.csr_matrix((self.hop['t'].conj(), (self.hop['j'], self.hop['i'])), 
                                                        shape=(self.lat.sites, self.lat.sites))

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
            self.pola = np.zeros((self.lat.sites, len(self.lat.tags)))
            for i, t in enumerate(self.lat.tags):
                self.pola[:, i] = np.sum(np.abs(self.vn[self.lat.coor['tag'] == t, :]) ** 2, axis=0)
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
        test_get_state_pola(self.vn, tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmax(self.pola[:, i_tag])
        print('State with polarization:', self.pola[ind, i_tag])
        return self.intensity[:, ind]

    def get_states_en(self, en_lims):
        '''
        Get, if any, the intensity of the sum of the states 
        between *en_lims[0]* and *en_lims[1]*.

        :param en_lims: List, en_lims[0] energy min, en_lims[1] energy max.

        :returns:
            * **intensity** -- Sum of the intensities between *e_min* and *e_max*.
        '''
        test_get_states_en(self.vn, en_lims)
        ind = np.where((self.en > en_lims[0]) & (self.en < en_lims[1]))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), en_lims[0], en_lims[1]))
        intensity = np.sum(np.abs(self.vn[:, ind]) ** 2, axis=1)
        return intensity




