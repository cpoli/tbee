from latticeTB import latticeTB
from eigTB import eigTB
from plotTB import plotTB

import matplotlib.pyplot as plt
import unittest
import numpy as np

def init():
    lat = latticeTB([[0, 0]], [b'a'], 1., 1.6)
    lat.get_lattice(nx=5, ny=5)
    eig = eigTB(lat)
    return eig
        
class TestEigTB(unittest.TestCase):
    '''
    Unittest class of **eigTB**.
    '''

    def test_errors_print_hop(self):
        # dict_hop not a dictionary with 0 <keys< n_max and values not dictionaries.
        eig = init()
        self.assertRaises(TypeError, eig.print_hop, 1j)
        self.assertRaises(ValueError, eig.print_hop, 0)
        self.assertRaises(ValueError, eig.print_hop, 100)

    def test_errors_set_ons(self):
        # on must be a list of numbers.
        eig = init()
        self.assertRaises(TypeError, eig.set_ons, 0)
        self.assertRaises(ValueError, eig.set_ons, [b'a'])

    def test_errors_set_hop_uni(self):
        # dict_hop not a dictionary with 0 <keys< n_max and values not numbers.
        eig = init()
        self.assertRaises(TypeError, eig.set_hop_uni, {0, 1})
        self.assertRaises(ValueError, eig.set_hop_uni, {'a': 1})
        self.assertRaises(ValueError, eig.set_hop_uni, {1: 'a'})
        self.assertRaises(ValueError, eig.set_hop_uni, {100: 1})

    def test_errors_set_hop(self):
        # dict_hop not a dictionary with 0 <keys< n_max and values not dictionaries.
        eig = init()
        self.assertRaises(TypeError, eig.set_hop, 0)
        self.assertRaises(ValueError, eig.set_hop, {'a': 1})
        self.assertRaises(TypeError, eig.set_hop, {1: 'a'})
        self.assertRaises(TypeError, eig.set_hop, {1: [1, 1]})
        self.assertRaises(ValueError, eig.set_hop, {1: {-1: 1}})
        self.assertRaises(ValueError, eig.set_hop, {1: {360: 1}})

    def test_set_hop_nearest(self):
        lat = latticeTB([[0, 0], [1, 0], [0, 1]], [b'a', b'b', b'c'], 2., np.pi/2)
        lat.get_lattice(nx=4, ny=4)
        eig = eigTB(lat)
        self.assertRaises(TypeError, eig.set_hop_nearest, 0)
        self.assertRaises(ValueError, eig.set_hop_nearest, {'ab': 2, b'ba': 1, b'ac': 2, b'ca': 1})
        self.assertRaises(ValueError, eig.set_hop_nearest, {b'ab': 'a', b'ba': 1, b'ac': 2, b'ca': 1})

    def test_set_dimerization_def(self):
        lat = latticeTB([[0, 0], [1, 0], [0, 1]], [b'a', b'b', b'c'], 2., np.pi/2)
        lat.get_lattice(nx=4, ny=4)
        eig = eigTB(lat)
        self.assertRaises(RuntimeError, eig.set_dimerization_def,{b'ab': 2, b'ba': 1, b'ac': 2, b'ca': 1}, 0, 0)
        eig.set_hop_nearest({b'ab': 2, b'ba': 1, b'ac': 2, b'ca': 1})
        self.assertRaises(TypeError, eig.set_dimerization_def, {b'ab': 'a', b'ba': 1, b'ac': 2, b'ca': 1}, 4, 0)
        self.assertRaises(TypeError, eig.set_dimerization_def, {b'ab': 2, b'ba': 1, b'ac': 2, b'ca': 1}, 'a', 0)
        self.assertRaises(TypeError, eig.set_dimerization_def, {b'ab': 2, b'ba': 1, b'ac': 2, b'ca': 1}, 4, 'a')
        self.assertRaises(ValueError, eig.set_dimerization_def, {b'asb': 2, b'ba': 1, b'ac': 2, b'ca': 1}, 4, 0)
        self.assertRaises(ValueError, eig.set_dimerization_def, {b'zb': 2, b'ba': 1, b'ac': 2, b'ca': 1}, 4, 0)

    def test_errors_set_hop_def(self):
        eig = init()
        self.assertRaises(RuntimeError, eig.set_hop_def, {(0, 1): 0})
        eig.set_hop_uni({1: 1})
        self.assertRaises(TypeError, eig.set_hop_def, 0)
        self.assertRaises(TypeError, eig.set_hop_def,{(1, 0): 'a'})
        self.assertRaises(ValueError, eig.set_hop_def, {(0, 0): 1})
        #self.assertRaises(ValueError, eig.set_hop_def, {(-1, 0): 1})
        #self.assertRaises(ValueError, eig.set_hop_def,{(100, 0): 1})


    def test_errors_set_ons_def(self):
        eig = init()
        self.assertRaises(TypeError, eig.set_ons_def, 0)
        self.assertRaises(ValueError, eig.set_ons_def, {-1: 0})
        self.assertRaises(TypeError, eig.set_ons_def, {0: 'a'})


    def test_errors_set_hop_disorder(self):
        eig = init()
        self.assertRaises(RuntimeError, eig.set_hop_disorder, 1)
        eig.set_hop_uni({1: 1})
        # alpha must be a positive number.
        self.assertRaises(TypeError, eig.set_hop_disorder, 0j)
        self.assertRaises(ValueError, eig.set_hop_disorder, 0)

    def test_errors_set_disorder_on(self):
        eig = init()
        # alpha must be a positive number.
        self.assertRaises(TypeError, eig.set_ons_disorder, 0j)
        self.assertRaises(ValueError, eig.set_ons_disorder, 0)

    def test_errors_get_ham(self):
        # complex_transpose must be a bool.
        eig = init()
        self.assertRaises(RuntimeError, eig.get_ham)
        eig.set_hop_uni({1: 1})
        self.assertRaises(TypeError, eig.get_ham, complex_transpose='a')

    def test_errors_get_eig(self):
        # eigenvec must be a bool.
        eig = init()
        eig.set_hop_uni({1: 1})
        self.assertRaises(RuntimeError, eig.get_eig)
        eig.get_ham()
        self.assertRaises(TypeError, eig.get_eig, eigenvec='a')

    def test_errors_get_state_pola(self):
        # tag_pola must be a binary char belonging to tags.
        eig = init()
        eig.set_hop_uni({1: 1})
        eig.get_ham()
        self.assertRaises(RuntimeError, eig.get_state_pola, tag_pola=b'a')
        eig.get_eig()
        self.assertRaises(RuntimeError, eig.get_state_pola, tag_pola=b'a')
        eig.get_eig(eigenvec=True)
        self.assertRaises(TypeError, eig.get_state_pola, tag_pola=0)
        self.assertRaises(ValueError, eig.get_state_pola, tag_pola=b'z')

    def test_errors_get_state_en(self):
        # tag_pola must be a binary char belonging to tags.
        eig = init()
        eig.set_hop_uni({1: 1})
        eig.get_ham()
        self.assertRaises(RuntimeError, eig.get_states_en, e_min=1., e_max='a')
        eig.get_eig(eigenvec=True)
        self.assertRaises(TypeError, eig.get_states_en, e_min='a', e_max=1.)
        self.assertRaises(TypeError, eig.get_states_en, e_min=1., e_max='a')
        self.assertRaises(ValueError, eig.get_states_en, e_min=2., e_max=1.)


if __name__ == '__main__':
    unittest.main()
