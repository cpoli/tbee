from tbee.lattice import Lattice
from tbee.system import System
from tbee.plot import Plot, plot
from tbee.graphene import GrapheneLattice, GrapheneSystem
import unittest
import numpy as np
import matplotlib.pyplot as plt


def build_system():
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (0.5, 0.)}]
    prim_vec = [(1., 0.), (0., 1.)]
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=4, n2=4)
    sys = System(lat)
    sys.set_hopping([{'n': 1, 't': 1.}])
    sys.set_hopping([{'n': 2, 't': 0.5}], upper_part=False)
    sys.set_onsite({'a': 0., 'b': 0.1})
    sys.get_ham()
    sys.get_eig(eigenvec=True, left=True)
    sys.get_ipr()
    sys.get_petermann()
    sys.get_coor_hop()
    return sys


def build_butterfly_system():
    lat = GrapheneLattice()
    lat.hexagon_zigzag(2)
    sys = GrapheneSystem(lat)
    sys.get_butterfly(t=1., N=5)
    return sys


class TestPlot(unittest.TestCase):
    '''
    Unittest of class **Plot**.
    '''
    def tearDown(self):
        plt.close('all')

    def test_init_error(self):
        self.assertRaises(TypeError, Plot, sys=0)

    def test_backward_compatible_alias(self):
        self.assertIs(plot, Plot)

    def test_lattice(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.lattice()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.lattice(plt_hop=True, plt_hop_low=True, plt_index=True, axis=True)
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_lattice_hop(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.lattice_hop()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.lattice_hop(plt_hop=True, plt_hop_low=True, plt_index=True)
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_lattice_custom_colors(self):
        sys = build_system()
        p = Plot(sys, colors=['orange', 'purple'])
        fig = p.lattice()
        self.assertEqual(fig.__class__.__name__, 'Figure')

    def test_spectrum_hist(self):
        sys = build_system()
        p = Plot(sys)
        p.spectrum_hist()
        p.spectrum_hist(lims=[-1., 1.])

    def test_spectrum_plain(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.spectrum()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.spectrum(lims=[-0.5, 0.5])
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_spectrum_tag_pola(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.spectrum(tag_pola='a')
        self.assertEqual(fig.__class__.__name__, 'Figure')

    def test_spectrum_ipr(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.spectrum(ipr=True)
        self.assertEqual(fig.__class__.__name__, 'Figure')

    def test_spectrum_peterman(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.spectrum(peterman=True)
        self.assertEqual(fig.__class__.__name__, 'Figure')

    def test_polarization_standalone(self):
        sys = build_system()
        p = Plot(sys)
        fig, ax = p.polarization(tag_pola='a')
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2, ax2 = p.polarization(tag_pola='a', lims=[-0.5, 0.5])
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_ipr_standalone(self):
        sys = build_system()
        p = Plot(sys)
        fig, ax = p.ipr()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2, ax2 = p.ipr(lims=[-0.5, 0.5])
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_petermann_standalone(self):
        sys = build_system()
        p = Plot(sys)
        fig, ax = p.petermann()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2, ax2 = p.petermann(lims=[-0.5, 0.5])
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_spectrum_complex(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.spectrum_complex()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.spectrum_complex(lims=[-1., 1.])
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_intensity_1d(self):
        sys = build_system()
        p = Plot(sys)
        intensity = np.abs(sys.rn[:, 0]) ** 2
        fig = p.intensity_1d(intensity)
        self.assertEqual(fig.__class__.__name__, 'Figure')

    def test_intensity_disk(self):
        sys = build_system()
        p = Plot(sys)
        intensity = np.abs(sys.rn[:, 0]) ** 2
        fig = p.intensity_disk(intensity)
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.intensity_disk(intensity, lims=[0., 1.], figsize=(4, 4))
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_intensity_area(self):
        sys = build_system()
        p = Plot(sys)
        intensity = np.abs(sys.rn[:, 0]) ** 2
        fig = p.intensity_area(intensity)
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.intensity_area(intensity, plt_hop=True)
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_butterfly(self):
        sys = build_butterfly_system()
        p = Plot(sys)
        fig = p.butterfly(sys.betas, sys.butterfly)
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.butterfly(sys.betas, sys.butterfly, lims=[-2., 2.], title='test')
        self.assertEqual(fig2.__class__.__name__, 'Figure')

    def test_show(self):
        sys = build_system()
        p = Plot(sys)
        p.show()

    def test_dos(self):
        sys = build_system()
        p = Plot(sys)
        fig = p.dos()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = p.dos(kernel='lorentzian', broadening=0.1,
                            e_grid=np.linspace(-3., 3., 101))
        self.assertEqual(fig2.__class__.__name__, 'Figure')


if __name__ == '__main__':
    unittest.main()
