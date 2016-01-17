from TB.lattice import *
from TB.system import *

import unittest
import numpy as np

def init():
    unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
    prim_vec = [(1., 0.), (0, 1.)]
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

    def test_print_distances(self):
        sys = init()
        self.assertRaises(TypeError, sys.print_distances, 1j)
        self.assertRaises(ValueError, sys.print_distances, 0)
        self.assertRaises(ValueError, sys.print_distances, 100)

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

    def test_set_hopping_example(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}, 
                          {'tag': b'b', 'r0': (1, 0)}]
        prim_vec = [(2, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        sys = system(lat)
        n1, n2 = 2, 2
        lat.get_lattice(n1=n1, n2=n2)
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 'ang': 10., 't': 1}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 4, 'tag': b'ab', 't': 1}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 'tag': b'ab', 'ang': 10., 't': 1}])
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 'ang': 90., 't': 1}], low=True)
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 'ang': 10., 't': 1}], low=True)
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 4, 'tag': b'ab', 't': 1}], low=True)
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 1, 'tag': b'ab', 'ang': 10., 't': 1}], low=True)
        self.assertRaises(ValueError, sys.set_hopping, [{'n': 4, 'tag': b'zz', 't': 1}])
        # check n
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.get_ham()
        ham1 = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 't': 1.}], low=True)
        sys.get_ham()
        ham2 = sys.ham
        self.assertTrue(np.allclose(ham1.toarray(), ham2.toarray()))
        # check tag
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 'tag': b'aa', 't': 1.}])
        sys.set_hopping([{'n': 1, 'tag': b'aa', 't': 2.}])
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': 3.}])
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': 4.}])
        sys.set_hopping([{'n': 1, 'tag': b'ba', 't': 5.}])
        sys.set_hopping([{'n': 1, 'tag': b'ba', 't': 6.}])
        sys.set_hopping([{'n': 1, 'tag': b'bb', 't': 7.}])
        sys.set_hopping([{'n': 1, 'tag': b'bb', 't': 8.}])
        sys.get_ham()
        ham1 = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 'tag': b'aa', 't': 1.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'aa', 't': 2.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'ba', 't': 3.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'ba', 't': 4.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': 5.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': 6.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'bb', 't': 7.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'bb', 't': 8.}], low=True)
        sys.get_ham()
        ham2 = sys.ham
        ham = [[ 0.+0.j,  4.+0.j,  0.+0.j,  0.+0.j,  2.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],
                    [ 4.+0.j,  0.+0.j,  6.+0.j,  0.+0.j,  0.+0.j,  8.+0.j,  0.+0.j,  0.+0.j],
                    [ 0.+0.j,  6.+0.j,  0.+0.j,  4.+0.j,  0.+0.j,  0.+0.j,  2.+0.j,  0.+0.j],
                    [ 0.+0.j,  0.+0.j,  4.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  8.+0.j],
                    [ 2.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  0.+0.j,  4.+0.j,  0.+0.j,  0.+0.j],
                    [ 0.+0.j,  8.+0.j,  0.+0.j,  0.+0.j,  4.+0.j,  0.+0.j,  6.+0.j,  0.+0.j],
                    [ 0.+0.j,  0.+0.j,  2.+0.j,  0.+0.j,  0.+0.j,  6.+0.j,  0.+0.j,  4.+0.j],
                    [ 0.+0.j,  0.+0.j,  0.+0.j,  8.+0.j,  0.+0.j,  0.+0.j,  4.+0.j,  0.+0.j]]
        self.assertTrue(np.allclose(ham, ham1.toarray()))
        self.assertTrue(np.allclose(ham, ham2.toarray()))
        # check angle
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 'ang': 0, 't': 1.}])
        sys.set_hopping([{'n': 1, 'ang': 90, 't': 2.}])
        sys.set_hopping([{'n': 1, 'ang': 0, 't': 3.}])
        sys.set_hopping([{'n': 1, 'ang': 90, 't': 4.}])
        sys.get_ham()
        ham1 = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 'ang': -180, 't': 1.}], low=True)
        sys.set_hopping([{'n': 1, 'ang': -90, 't': 2.}], low=True)
        sys.set_hopping([{'n': 1, 'ang': -180, 't': 3.}], low=True)
        sys.set_hopping([{'n': 1, 'ang': -90, 't': 4.}], low=True)
        sys.get_ham()
        ham2 = sys.ham
        self.assertTrue(np.allclose(ham1.toarray(), ham2.toarray()))
        sys.clear_hopping()
        sys.set_hopping([{'n': 4, 'ang': 26.566, 't': 1.}])
        sys.set_hopping([{'n': 4, 'ang': 26.564, 't': 2.}])
        sys.set_hopping([{'n': 4, 'ang': 153.433, 't': 3.}])
        sys.set_hopping([{'n': 4, 'ang': 153.435, 't': 4.}])
        sys.get_ham()
        ham1 = sys.ham
        sys.set_hopping([{'n': 4, 'ang': 26.565-180, 't': 2.}], low=True)
        sys.set_hopping([{'n': 4, 'ang': 153.434-180, 't': 4.}], low=True)
        sys.get_ham()
        ham2 = sys.ham
        self.assertTrue(np.allclose(ham1.toarray(), ham2.toarray()))
        # check tag & angle
        sys.clear_hopping()
        sys.set_hopping([{'n': 2, 'ang': 135, 'tag':  b'ab', 't': 1.}])
        sys.set_hopping([{'n': 2, 'ang': 135, 'tag':  b'ab', 't': 2.}])
        sys.set_hopping([{'n': 2, 'ang': 45, 'tag':  b'ab', 't': 2.}])
        sys.set_hopping([{'n': 2, 'ang': 45, 'tag':  b'ab', 't': 3.}])
        sys.set_hopping([{'n': 2, 'ang': 135, 'tag':  b'ba', 't': 3.}])
        sys.set_hopping([{'n': 2, 'ang': 135, 'tag':  b'ba', 't': 4.}])
        sys.set_hopping([{'n': 2, 'ang': 45, 'tag':  b'ba', 't': 4.}])
        sys.set_hopping([{'n': 2, 'ang': 45, 'tag':  b'ba', 't': 5.}])
        sys.get_ham()
        ham1 = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 2, 'ang': -45, 'tag':  b'ba', 't': 1.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -45, 'tag':  b'ba', 't': 2.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -135, 'tag':  b'ba', 't': 2.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -135, 'tag':  b'ba', 't': 3.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -45, 'tag':  b'ab', 't': 3.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -45, 'tag':  b'ab', 't': 4.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -135, 'tag':  b'ab', 't': 4.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -135, 'tag':  b'ab', 't': 5.}], low=True)
        sys.get_ham()
        ham2 = sys.ham
        self.assertTrue(np.allclose(ham1.toarray(), ham2.toarray()))


    def test_set_hopping_left(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}, 
                           {'tag': b'b', 'r0': (1, 0)}]
        prim_vec = [(2, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        sys = system(lat)
        n1, n2 = 2, 2
        lat.get_lattice(n1=n1, n2=n2)
        no_n1 = 10
        no_n2 = 6
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.set_hopping([{'n': 1, 't': 2.}])
        sys.set_hopping([{'n': 1, 't': 3.}])
        sys.set_hopping([{'n': 1, 't': 4.}])
        sys.set_hopping([{'n': 1, 'tag': b'ab', 't': 5.}])
        self.assertTrue(len(sys.hop) == 10)
        sys.set_hopping([{'n': 1, 't': 4.}], low=True)
        sys.set_hopping([{'n': 1, 'tag': b'ba', 't': 5.}], low=True)
        self.assertTrue(len(sys.hop) == 20)
        self.assertTrue(np.all(sys.hop['i'][:10] == sys.hop['j'][10:]))
        self.assertTrue(np.all(sys.hop['j'][:10] == sys.hop['i'][10:]))
        self.assertTrue(np.all(sys.hop['ang'][:10] == sys.hop['ang'][10:]+180))
        self.assertTrue(np.all(sys.hop['t'][:10] == sys.hop['t'][10:]))
        sys.set_hopping([{'n': 2, 'ang': 45, 't': 8.}])
        sys.set_hopping([{'n': 2, 'ang': 135, 't': 0.}])
        sys.set_hopping([{'n': 2, 'ang': 135, 't': 9.}])
        self.assertTrue(len(sys.hop) == 26)
        sys.set_hopping([{'n': 2, 'ang': -135, 't': 8.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -45, 't': 2.}], low=True)
        sys.set_hopping([{'n': 2, 'ang': -45, 't': 9.}], low=True)
        self.assertTrue(len(sys.hop) == 32)
        self.assertTrue(np.all(sys.hop['i'][20: 26] == sys.hop['j'][26:]))
        self.assertTrue(np.all(sys.hop['j'][20: 26] == sys.hop['i'][26:]))
        self.assertTrue(np.all(sys.hop['ang'][20: 26] == sys.hop['ang'][26:]+180))
        self.assertTrue(np.all(sys.hop['t'][20: 26] == sys.hop['t'][26:]))
        sys.set_hopping([{'n': 4, 'ang': 26.565, 'tag': b'aa', 't': 11}])
        sys.set_hopping([{'n': 4, 'ang': 26.565, 'tag': b'aa', 't': 99}])
        self.assertTrue(len(sys.hop) == 33)
        sys.set_hopping([{'n': 4, 'ang': 26.565-180, 'tag': b'aa', 't': 22}], low=True)
        sys.set_hopping([{'n': 4, 'ang': 26.565-180, 'tag': b'aa', 't': 99}], low=True)
        self.assertTrue(len(sys.hop) == 34)
        self.assertTrue(sys.hop['i'][32] == sys.hop['j'][33])
        self.assertTrue(sys.hop['j'][33] == sys.hop['i'][32])
        self.assertTrue(np.isclose(sys.hop['ang'][32], sys.hop['ang'][33]+180.))
        self.assertTrue(sys.hop['t'][32] == sys.hop['t'][33])

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
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}, {'tag': b'b', 'r0': (1, 0)}]
        prim_vec = [(2., 0.)]
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
        self.assertTrue(np.allclose(sys.en, en) == True)
        self.assertTrue(np.allclose(sys.get_intensity_pola_max(tag_pola=b'a'), 
                                                   zero_mode) == True)
        self.assertTrue(np.allclose(sys.get_intensity_en(lims=[-3, 1e-3]), 
                                                   modes_neg) == True)

    def test_square(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1., 0.), (0., 1.)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1, n2 = 3, 3
        lat.get_lattice(n1=n1, n2=n2)
        sys = system(lat=lat)
        t_1, t_2, t_3, t_4, t_5 = 1., 2., 3., 4., 5. 
        sys.set_hopping([{'n': 1, 't': t_1}, 
                                    {'n': 2, 't': t_2}, 
                                    {'n': 3, 't': t_3},
                                    {'n': 4, 't': t_4},
                                    {'n': 5, 't': t_5}])
        self.assertTrue(np.sum(sys.hop['n'] == 1) == 12)
        self.assertTrue(np.sum(sys.hop['n'] == 2) == 8)
        self.assertTrue(np.sum(sys.hop['n'] == 3) == 6)
        self.assertTrue(np.sum(sys.hop['n'] == 4) == 8)
        self.assertTrue(np.sum(sys.hop['n'] == 5) == 2)
        sys.get_ham()
        ham_up = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 't': t_1}, 
                                    {'n': 2, 't': t_2}, 
                                    {'n': 3, 't': t_3},
                                    {'n': 4, 't': t_4},
                                    {'n': 5, 't': t_5}], low=True)
        self.assertTrue(np.sum(sys.hop['n'] == 1) == 12)
        self.assertTrue(np.sum(sys.hop['n'] == 2) == 8)
        self.assertTrue(np.sum(sys.hop['n'] == 3) == 6)
        self.assertTrue(np.sum(sys.hop['n'] == 4) == 8)
        self.assertTrue(np.sum(sys.hop['n'] == 5) == 2)
        sys.get_ham()
        ham_low = sys.ham
        sys.clear_hopping()
        sys.set_hopping([{'n': 1, 't': t_1}, 
                                    {'n': 2, 't': t_2}, 
                                    {'n': 3, 't': t_3},
                                    {'n': 4, 't': t_4},
                                    {'n': 5, 't': t_5}])        
        sys.set_hopping([{'n': 1, 't': t_1}, 
                                    {'n': 2, 't': t_2}, 
                                    {'n': 3, 't': t_3},
                                    {'n': 4, 't': t_4},
                                    {'n': 5, 't': t_5}], low=True)
        self.assertTrue(np.sum(sys.hop['n'] == 1) == 24)
        self.assertTrue(np.sum(sys.hop['n'] == 2) == 16)
        self.assertTrue(np.sum(sys.hop['n'] == 3) == 12)
        self.assertTrue(np.sum(sys.hop['n'] == 4) == 16)
        self.assertTrue(np.sum(sys.hop['n'] == 5) == 4)
        sys.get_ham()
        ham_full = sys.ham
        ham = [ [ 0.,  1.,  3.,  1.,  2.,  4.,  3.,  4.,  5.],
                     [ 1.,  0.,  1.,  2.,  1.,  2.,  4.,  3.,  4.],
                     [ 3.,  1.,  0.,  4.,  2.,  1.,  5.,  4.,  3.],
                     [ 1.,  2.,  4.,  0.,  1.,  3.,  1.,  2.,  4.],
                     [ 2.,  1.,  2.,  1.,  0.,  1.,  2.,  1.,  2.],
                     [ 4.,  2.,  1.,  3.,  1.,  0.,  4.,  2.,  1.],
                     [ 3.,  4.,  5.,  1.,  2.,  4.,  0.,  1.,  3.],
                     [ 4.,  3.,  4.,  2.,  1.,  2.,  1.,  0.,  1.],
                     [ 5.,  4.,  3.,  4.,  2.,  1.,  3.,  1.,  0.], ]
        self.assertTrue(np.allclose(ham, ham_up.toarray()))
        self.assertTrue(np.allclose(ham, ham_low.toarray()))
        self.assertTrue(np.allclose(ham, ham_full.toarray()))

if __name__ == '__main__':
    unittest.main()
