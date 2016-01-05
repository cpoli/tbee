from lattice import *
from system import *

import unittest
import numpy as np

def init():
    unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
    prim_vec = {'norm': 1, 'angle': 0}
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=5, n2=5)
    sys = system(lat)
    return sys
        
class TestSystem(unittest.TestCase):
    '''
    Unittest of class **system**.
    '''
    def test_init(self):
        self.assertRaises(TypeError, system, lat=0)

    def test_print_hopping(self):
        sys = init()
        self.assertRaises(TypeError, sys.print_hopping, 1j)
        self.assertRaises(ValueError, sys.print_hopping, 0)
        self.assertRaises(ValueError, sys.print_hopping, 100)

    def test_set_onsite(self):
        sys = init()
        self.assertRaises(TypeError, sys.set_onsite, 0)
        self.assertRaises(ValueError, sys.set_onsite, {b'a': 'a'})
        self.assertRaises(ValueError, sys.set_onsite, {b'z': 1.})

    def test_set_hopping(self):
        sys = init()
        self.assertRaises(TypeError, sys.set_hopping, 0)
        self.assertRaises(KeyError, sys.set_hopping, [{'n': 1}])
        # keys "n" and "t"
        self.assertRaises(TypeError, sys.set_hopping, [{'n': 'a', 't': 1.}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': -1, 't': 1.}])
        self.assertRaises(TypeError, sys.set_hopping, [{'n': 1, 't': 'a'}])
        # keys "ang"
        self.assertRaises(TypeError, sys.set_hopping, [{'n': 1, 't': 1., 'ang': 'a'}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 't': 1., 'ang': 360}])
        # keys "tag"
        self.assertRaises(TypeError, sys.set_hopping, [{'n': 1, 't': 1., 'tag': 0}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 't': 1., 'tag': b'a'}])

        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 't': 1., 'tag': b'aa', 'ang': 0, 'a':0}])

    def test_get_ham(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.get_ham)

    def test_get_eig(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.get_eig)
        sys.set_hopping([{'n': 1, 't': 1.}])
        self.assertRaises(RuntimeError, sys.get_eig)
        sys.get_ham()
        self.assertRaises(TypeError, sys.get_eig, eigenvec='a')

    def test_get_intensity_pola_max(self):
        sys = init()
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.get_ham()
        self.assertRaises(RuntimeError, sys.get_intensity_pola_max, tag_pola=b'a')
        sys.get_eig()
        self.assertRaises(RuntimeError, sys.get_intensity_pola_max, tag_pola=b'a')
        sys.get_eig(eigenvec=True)
        self.assertRaises(TypeError, sys.get_intensity_pola_max, tag_pola=0)
        self.assertRaises(ValueError, sys.get_intensity_pola_max, tag_pola=b'z')

    def test_get_intensity_pola_min(self):
        sys = init()
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.get_ham()
        self.assertRaises(RuntimeError, sys.get_intensity_pola_max, tag_pola=b'b')
        sys.get_eig()
        self.assertRaises(RuntimeError, sys.get_intensity_pola_max, tag_pola=b'b')
        sys.get_eig(eigenvec=True)
        self.assertRaises(TypeError, sys.get_intensity_pola_max, tag_pola=0)
        self.assertRaises(ValueError, sys.get_intensity_pola_max, tag_pola=b'b')

    def test_set_hopping_def(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.set_hopping_def, {(0, 1): 0})
        sys.set_hopping([{'n':1, 't':1.}])
        self.assertRaises(TypeError, sys.set_hopping_def, 0)
        self.assertRaises(TypeError, sys.set_hopping_def, {(1, 0): 'a'})
        self.assertRaises(ValueError, sys.set_hopping_def, {(0, 0): 1})

    def test_set_onsite_def(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.set_onsite_dis, 0.1)
        sys.set_onsite({b'a': 0.})
        self.assertRaises(TypeError, sys.set_onsite_def, 0)
        self.assertRaises(TypeError, sys.set_onsite_def, {0: 'a'})

    def test_set_onsite_dis(self):
        sys = init()
        # alpha must be a positive number.
        self.assertRaises(RuntimeError, sys.set_onsite_dis, 0.1)
        sys.set_onsite({b'a': 0.})
        self.assertRaises(TypeError, sys.set_onsite_dis, 'a')

    def test_set_hopping_dis(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.set_hopping_dis, 1.)
        sys.set_hopping([{'n': 1, 't': 1.}])
        # alpha must be a positive number.
        self.assertRaises(TypeError, sys.set_hopping_dis, 'a')

    def test_get_coor_hop(self):
        sys = init()
        self.assertRaises(RuntimeError, sys.get_coor_hop)
        sys.set_hopping([{'n': 2, 't': 1.}])
        self.assertRaises(ValueError, sys.get_coor_hop)

    def test_get_intensity_en(self):
        # tag_pola must be a binary char belonging to lat.tags.
        sys = init()
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.get_ham()
        self.assertRaises(RuntimeError, sys.get_intensity_en, lims=[1., 2.])
        sys.get_eig(eigenvec=True)
        self.assertRaises(TypeError, sys.get_intensity_en, lims=['a', 2.])
        self.assertRaises(TypeError, sys.get_intensity_en, lims=[1., 'a'])
        self.assertRaises(ValueError, sys.get_intensity_en, lims=[2., 1.])

    def test_dimer_chain(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}, {'tag': b'b', 'r0': [1, 0]}]
        prim_vec = {'norm': 2, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1 = 10
        lat.get_lattice(n1=n1)
        lat.remove_sites(index=[2*n1-1])
        sys = system(lat=lat)
        e_a, e_b = 0j, -1j 
        t_ab, t_ba = 2., 1. 
        sys.set_onsite({b'a': e_a, b'b': e_b})
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': t_ab}, 
                                    {'n': 1, 'tag': b'ba', 't': t_ba}])
        sys.get_ham()
        sys.get_eig(eigenvec=True)
        onsite = [ 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0., -1j, 0]
        hop = np.array([(1, 0, 1, (2+0j), 0, b'ab'), (1, 1, 2, (1+0j), 0, b'ba'),
                                 (1, 2, 3, (2+0j), 0, b'ab'), (1, 3, 4, (1+0j), 0, b'ba'),
                                 (1, 4, 5, (2+0j), 0, b'ab'), (1, 5, 6, (1+0j), 0, b'ba'),
                                 (1, 6, 7, (2+0j), 0, b'ab'), (1, 7, 8, (1+0j), 0, b'ba'),
                                 (1, 8, 9, (2+0j), 0, b'ab'), (1, 9, 10, (1+0j), 0, b'ba'),
                                 (1, 10, 11, (2+0j), 0, b'ab'), (1, 11, 12, (1+0j), 0, b'ba'),
                                 (1, 12, 13, (2+0j), 0, b'ab'), (1, 13, 14, (1+0j), 0, b'ba'),
                                 (1, 14, 15, (2+0j), 0, b'ab'), (1, 15, 16, (1+0j), 0, b'ba'),
                                 (1, 16, 17, (2+0j), 0, b'ab'), (1, 17, 18, (1+0j), 0, b'ba')], 
                                 dtype=[('n', 'u2'), ('i', 'u4'), ('j', 'u4'), 
                                             ('t', 'c16'), ('ang', 'i2'), ('tag', 'S2')])
        en = [ -2.92476086e+00 -5.00000000e-01j,
                  -2.82596319e+00 -5.00000000e-01j,
                  -2.66479662e+00 -5.00000000e-01j,
                  -2.44664423e+00 -5.00000000e-01j,
                  -2.17944947e+00 -5.00000000e-01j,
                  -1.87454848e+00 -5.00000000e-01j,
                  -1.54882504e+00 -5.00000000e-01j,
                  -1.23041945e+00 -5.00000000e-01j,
                  -9.72509092e-01 -5.00000000e-01j,
                  -8.62722135e-16 -1.59417370e-16j,
                   9.72509092e-01 -5.00000000e-01j,
                   1.23041945e+00 -5.00000000e-01j,
                   1.54882504e+00 -5.00000000e-01j,
                   1.87454848e+00 -5.00000000e-01j,
                   2.17944947e+00 -5.00000000e-01j,
                   2.44664423e+00 -5.00000000e-01j,
                   2.66479662e+00 -5.00000000e-01j,
                   2.82596319e+00 -5.00000000e-01j,
                   2.92476086e+00 -5.00000000e-01j]
        zero_mode = [  2.86102568e-06,   1.27709173e-32,   
                                 1.14441027e-05,   3.16009683e-32,   
                                 4.57764108e-05,   7.31793779e-33,
                                 1.83105643e-04,   2.63340172e-32,   
                                 7.32422573e-04,   3.41224428e-32,   
                                 2.92969029e-03,   1.73027735e-33,
                                 1.17187612e-02,   1.14798387e-32,
                                 4.68750447e-02,   5.25130197e-34,   
                                 1.87500179e-01,   2.41658191e-33,
                                 7.50000715e-01]
        modes_neg = [ 0.50000143,  0.5      ,  
                                 0.50000572,  0.5      ,  
                                0.50002289,   0.5      ,  
                                0.50009155,  0.5       ,  
                                0.50036621,  0.5       ,
                                0.50146485,  0.5       ,  
                                0.50585938,  0.5       ,  
                                0.52343752,  0.5       ,  
                                0.59375009,  0.5       ,  
                                0.87500036]
        self.assertTrue(np.allclose(sys.onsite, onsite) == True)
        self.assertTrue(np.allclose(sys.hop['i'], hop['i']) == True)
        self.assertTrue(np.allclose(sys.hop['j'], hop['j']) == True)
        self.assertTrue(np.allclose(sys.hop['t'], hop['t']) == True)
        self.assertTrue(np.array_equal(sys.hop['tag'], hop['tag']) == True)
        self.assertTrue(np.allclose(sys.en, en) == True)
        self.assertTrue(np.allclose(sys.get_intensity_pola_max(tag_pola=b'a'), 
                                                   zero_mode) == True)
        self.assertTrue(np.allclose(sys.get_intensity_en(lims=[-3, 1e-3]), 
                                                   modes_neg) == True)

if __name__ == '__main__':
    unittest.main()
