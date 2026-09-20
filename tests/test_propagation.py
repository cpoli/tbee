from tbee.lattice import lattice
from tbee.system import system
from tbee.propagation import Propagation, propagation
import unittest
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from unittest import mock


PI = np.pi


def build_chain():
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)},
                      {'tag': 'b', 'r0': (0.5, 0.5/sqrt(3))}]
    prim_vec = [(1, 0.),
                (np.cos(PI/3), np.sin(PI/3))]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    n1, n2 = 2, 2
    lat.get_lattice(n1=n1, n2=n2)
    lat.remove_dangling()
    sys = system(lat=lat)
    sys.set_hopping([{'n': 1, 't': 1.}])
    sys.get_ham()
    return sys


def build_propagated(norm=True, steps=20):
    sys = build_chain()
    prop = Propagation(sys.lat)
    psi_init = np.zeros(sys.lat.sites, 'c16')
    psi_init[0] = 1.
    prop.get_propagation(sys.ham, psi_init, steps=steps, dz=0.1, norm=norm)
    return prop


def build_dimer_propagated():
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (1., 0.)}]
    prim_vec = [(2., 0.)]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=1)
    sys = system(lat=lat)
    sys.set_hopping_manual({(0, 1): 1.})
    sys.get_ham()
    prop = Propagation(lat)
    psi_init = np.array([1., 0.], 'c16')
    prop.get_propagation(sys.ham, psi_init, steps=20, dz=0.1, norm=True)
    return prop


class TestPropagation(unittest.TestCase):
    '''
    Unittest of class **Propagation**.
    '''
    def tearDown(self):
        plt.close('all')

    def test_init_error(self):
        self.assertRaises(TypeError, Propagation, lat=0)

    def test_backward_compatible_alias(self):
        self.assertIs(propagation, Propagation)

    def test_get_propagation(self):
        prop = build_propagated(norm=True)
        self.assertEqual(prop.prop.shape, (prop.lat.sites, 20))
        prop2 = build_propagated(norm=False)
        self.assertEqual(prop2.prop.shape, (prop2.lat.sites, 20))

    def test_get_pumping(self):
        sys = build_chain()
        sys2 = build_chain()
        sys2.set_hopping([{'n': 1, 't': 2.}])
        sys2.get_ham()
        prop = Propagation(sys.lat)
        psi_init = np.zeros(sys.lat.sites, 'c16')
        psi_init[0] = 1.
        prop.get_pumping([sys.ham, sys2.ham], psi_init, steps=21, dz=0.1, norm=True)
        self.assertEqual(prop.prop.shape, (sys.lat.sites, 21))
        prop.get_pumping([sys.ham, sys2.ham], psi_init, steps=21, dz=0.1, norm=False)
        self.assertEqual(prop.prop.shape, (sys.lat.sites, 21))

    def test_plt_propagation_1d(self):
        prop = build_propagated()
        for prop_type in ('real', 'imag', 'norm'):
            fig = prop.plt_propagation_1d(prop_type=prop_type)
            self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = prop.plt_propagation_1d(figsize=(6, 4))
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_get_animation_posix(self):
        prop = build_propagated()
        for prop_type in ('real', 'imag', 'norm'):
            ani = prop.get_animation(prop_type=prop_type)
            self.assertEqual(ani.__class__.__name__, 'FuncAnimation')

    def test_get_animation_non_posix(self):
        prop = build_propagated()
        with mock.patch('os.name', 'nt'):
            ani = prop.get_animation()
        self.assertEqual(ani.__class__.__name__, 'FuncAnimation')

    def test_get_animation_nb(self):
        prop = build_propagated()
        ani = prop.get_animation_nb(prop_type='real')
        self.assertEqual(ani.__class__.__name__, 'FuncAnimation')
        # Force the (lazily-evaluated) per-frame callback to actually run.
        ani.to_jshtml()
        ani2 = prop.get_animation_nb(prop_type='norm')
        self.assertEqual(ani2.__class__.__name__, 'FuncAnimation')

    def test_plt_prop_dimer(self):
        prop = Propagation(build_chain().lat)
        self.assertRaises(Exception, prop.plt_prop_dimer)
        prop_dimer = build_dimer_propagated()
        fig = prop_dimer.plt_prop_dimer()
        self.assertEqual(fig.__class__.__name__, 'Figure')


if __name__ == '__main__':
    unittest.main()
