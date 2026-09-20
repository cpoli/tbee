import unittest
import numpy as np
import scipy.sparse as sparse

import tbee.error_handling as eh
from tbee.lattice import Lattice, COOR_DTYPE
from tbee.system import System


class TestGeneric(unittest.TestCase):
    def test_boolean(self):
        eh.boolean(True, 'x')
        self.assertRaises(TypeError, eh.boolean, 1, 'x')

    def test_positive_int(self):
        eh.positive_int(1, 'x')
        self.assertRaises(TypeError, eh.positive_int, 1.0, 'x')
        self.assertRaises(ValueError, eh.positive_int, 0, 'x')

    def test_positive_int_lim(self):
        eh.positive_int_lim(1, 'x', 5)
        self.assertRaises(TypeError, eh.positive_int_lim, 1.0, 'x', 5)
        self.assertRaises(ValueError, eh.positive_int_lim, 0, 'x', 5)
        self.assertRaises(ValueError, eh.positive_int_lim, 6, 'x', 5)

    def test_real_number(self):
        eh.real_number(1.0, 'x')
        self.assertRaises(TypeError, eh.real_number, 'a', 'x')

    def test_positive_real(self):
        eh.positive_real(1.0, 'x')
        self.assertRaises(TypeError, eh.positive_real, 'a', 'x')
        self.assertRaises(ValueError, eh.positive_real, 0, 'x')

    def test_positive_real_zero(self):
        eh.positive_real_zero(0, 'x')
        self.assertRaises(TypeError, eh.positive_real_zero, 'a', 'x')
        self.assertRaises(ValueError, eh.positive_real_zero, -1, 'x')

    def test_negative_real(self):
        eh.negative_real(-1, 'x')
        self.assertRaises(TypeError, eh.negative_real, 'a', 'x')
        self.assertRaises(ValueError, eh.negative_real, 0, 'x')

    def test_number(self):
        eh.number(1j, 'x')
        self.assertRaises(TypeError, eh.number, 'a', 'x')

    def test_is_callable(self):
        eh.is_callable(lambda: None, 'x')
        self.assertRaises(TypeError, eh.is_callable, 0, 'x')

    def test_larger(self):
        eh.larger(1, 'a', 2, 'b')
        self.assertRaises(ValueError, eh.larger, 2, 'a', 1, 'b')

    def test_smaller(self):
        eh.smaller(1, 'a', 2, 'b')
        self.assertRaises(ValueError, eh.smaller, 2, 'a', 1, 'b')

    def test_string(self):
        eh.string(None, 'x')
        eh.string('a', 'x')
        self.assertRaises(TypeError, eh.string, 0, 'x')

    def test_ndarray(self):
        eh.ndarray(np.array([1, 2]), 'x', 2)
        self.assertRaises(TypeError, eh.ndarray, [1, 2], 'x', 2)
        self.assertRaises(ValueError, eh.ndarray, np.array([1, 2]), 'x', 3)

    def test_ndarray_null(self):
        eh.ndarray_null(np.array([1., 0.]), 'x')
        self.assertRaises(ValueError, eh.ndarray_null, np.array([0., 0.]), 'x')

    def test_ndarray_empty(self):
        eh.ndarray_empty(np.array([1.]), 'x')
        self.assertRaises(ValueError, eh.ndarray_empty, np.array([]), 'x')

    def test_list_tuple_2elem(self):
        eh.list_tuple_2elem(None, 'x')
        eh.list_tuple_2elem((1, 2), 'x')
        eh.list_tuple_2elem([1, 2], 'x')
        self.assertRaises(TypeError, eh.list_tuple_2elem, 'ab', 'x')
        self.assertRaises(ValueError, eh.list_tuple_2elem, (1, 2, 3), 'x')

    def test_tuple_2elem(self):
        eh.tuple_2elem(None, 'x')
        eh.tuple_2elem((1, 2), 'x')
        self.assertRaises(TypeError, eh.tuple_2elem, [1, 2], 'x')
        self.assertRaises(ValueError, eh.tuple_2elem, (1, 2, 3), 'x')


class TestLatticeChecks(unittest.TestCase):
    def test_lat(self):
        lat = Lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        eh.lat(lat)
        self.assertRaises(TypeError, eh.lat, 0)

    def test_unit_cell(self):
        eh.unit_cell([{'tag': 'a', 'r0': (0., 0.)}])
        self.assertRaises(TypeError, eh.unit_cell, 0)
        self.assertRaises(KeyError, eh.unit_cell, [{'r0': (0., 0.)}])
        self.assertRaises(KeyError, eh.unit_cell, [{'tag': 'a'}])
        self.assertRaises(TypeError, eh.unit_cell, [{'tag': 0, 'r0': (0., 0.)}])
        self.assertRaises(ValueError, eh.unit_cell, [{'tag': 'ab', 'r0': (0., 0.)}])
        self.assertRaises(TypeError, eh.unit_cell, [{'tag': 'a', 'r0': [0., 0.]}])
        self.assertRaises(ValueError, eh.unit_cell, [{'tag': 'a', 'r0': (0., 0., 0.)}])
        self.assertRaises(ValueError, eh.unit_cell, [{'tag': 'a', 'r0': (0., 'z')}])

    def test_prim_vec(self):
        eh.prim_vec([(1., 0.)])
        eh.prim_vec([(1., 0.), (0., 1.)])
        self.assertRaises(TypeError, eh.prim_vec, (1., 0.))
        self.assertRaises(ValueError, eh.prim_vec, [(1., 0.), (0., 1.), (1., 1.)])
        self.assertRaises(TypeError, eh.prim_vec, [[1., 0.]])
        self.assertRaises(ValueError, eh.prim_vec, [(1., 0., 0.)])
        self.assertRaises(ValueError, eh.prim_vec, [(1., 'z')])
        self.assertRaises(ValueError, eh.prim_vec, [(0.01, 0.01)])

    def test_get_lattice(self):
        eh.get_lattice([(1., 0.)], 1, 1)
        self.assertRaises(TypeError, eh.get_lattice, [(1., 0.)], 'a', 1)
        self.assertRaises(ValueError, eh.get_lattice, [(1., 0.)], 0, 1)
        self.assertRaises(ValueError, eh.get_lattice, [(1., 0.)], 1, 2)

    def test_coor(self):
        good = np.array([(0., 0., 'a')], dtype=COOR_DTYPE)
        eh.coor(good)
        bad = np.array([(0., 0.)], dtype=[('x', 'f8'), ('y', 'f8')])
        self.assertRaises(TypeError, eh.coor, bad)

    def test_empty_coor(self):
        eh.empty_coor(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_coor, np.array([]))

    def test_empty_coor_hop(self):
        eh.empty_coor_hop(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_coor_hop, np.array([]))

    def test_coor_1d(self):
        coor_flat = np.array([(0., 0., 'a'), (1., 0., 'a')], dtype=COOR_DTYPE)
        eh.coor_1d(coor_flat)
        coor_not_flat = np.array([(0., 0., 'a'), (1., 1., 'a')], dtype=COOR_DTYPE)
        self.assertRaises(ValueError, eh.coor_1d, coor_not_flat)

    def test_remove_sites(self):
        eh.remove_sites([0, 1], 5)
        self.assertRaises(TypeError, eh.remove_sites, 0, 5)
        self.assertRaises(ValueError, eh.remove_sites, [0.5], 5)
        self.assertRaises(ValueError, eh.remove_sites, [-1], 5)
        self.assertRaises(ValueError, eh.remove_sites, [5], 5)

    def test_shift(self):
        eh.shift(1.)
        self.assertRaises(TypeError, eh.shift, 'a')

    def test_boundary_line(self):
        eh.boundary_line(1., 1., 1.)
        self.assertRaises(TypeError, eh.boundary_line, 'a', 1., 1.)
        self.assertRaises(TypeError, eh.boundary_line, 1., 'a', 1.)
        self.assertRaises(TypeError, eh.boundary_line, 1., 1., 'a')

    def test_ellipse(self):
        eh.ellipse(1., 1.)
        self.assertRaises(TypeError, eh.ellipse, 'a', 1.)
        self.assertRaises(TypeError, eh.ellipse, 1., 'a')
        self.assertRaises(ValueError, eh.ellipse, -1., 1.)
        self.assertRaises(ValueError, eh.ellipse, 1., -1.)

    def test_sites(self):
        eh.sites(1)
        self.assertRaises(RuntimeError, eh.sites, 0)


class TestSystemChecks(unittest.TestCase):
    def test_sys(self):
        lat = Lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        lat.get_lattice(n1=2)
        sys = System(lat)
        eh.sys(sys)
        self.assertRaises(TypeError, eh.sys, 0)

    def test_print_hopping(self):
        eh.print_hopping(1, 5)
        self.assertRaises(TypeError, eh.print_hopping, 'a', 5)
        self.assertRaises(ValueError, eh.print_hopping, 0, 5)
        self.assertRaises(ValueError, eh.print_hopping, 5, 5)

    def test_set_onsite(self):
        eh.set_onsite({'a': 1.}, ['a', 'b'])
        self.assertRaises(TypeError, eh.set_onsite, 0, ['a'])
        self.assertRaises(ValueError, eh.set_onsite, {'z': 1.}, ['a'])
        self.assertRaises(ValueError, eh.set_onsite, {'a': 'z'}, ['a'])

    def test_set_hopping(self):
        eh.set_hopping([{'n': 1, 't': 1.}], 2)
        eh.set_hopping([{'n': 1, 't': 1., 'ang': 0.}], 2)
        eh.set_hopping([{'n': 1, 't': 1., 'tag': 'ab'}], 2)
        eh.set_hopping([{'n': 1, 't': 1., 'ang': 0., 'tag': 'ab'}], 2)
        self.assertRaises(TypeError, eh.set_hopping, 0, 2)
        self.assertRaises(TypeError, eh.set_hopping, [0], 2)
        self.assertRaises(KeyError, eh.set_hopping, [{'n': 1}], 2)
        self.assertRaises(TypeError, eh.set_hopping, [{'n': 'a', 't': 1.}], 2)
        self.assertRaises(ValueError, eh.set_hopping, [{'n': 0, 't': 1.}], 2)
        self.assertRaises(TypeError, eh.set_hopping, [{'n': 1, 't': 'a'}], 2)
        # len(dic) == 3, neither 'tag' nor 'ang': now correctly rejected.
        self.assertRaises(KeyError, eh.set_hopping, [{'n': 1, 't': 1., 'x': 0}], 2)
        # len(dic) == 4, neither 'tag' nor 'ang'.
        self.assertRaises(KeyError, eh.set_hopping, [{'n': 1, 't': 1., 'x': 0, 'y': 0}], 2)
        self.assertRaises(ValueError, eh.set_hopping,
                                    [{'n': 1, 't': 1., 'x': 0, 'y': 0, 'z': 0}], 2)
        self.assertRaises(TypeError, eh.set_hopping, [{'n': 1, 't': 1., 'tag': 0}], 2)
        self.assertRaises(ValueError, eh.set_hopping, [{'n': 1, 't': 1., 'tag': 'abc'}], 2)
        self.assertRaises(TypeError, eh.set_hopping, [{'n': 1, 't': 1., 'ang': 'a'}], 2)

    def test_index(self):
        eh.index(np.array([True, False]), {})
        self.assertRaises(ValueError, eh.index, np.array([False, False]), {})

    def test_set_hopping_def(self):
        eh.set_hopping_def(None, {(0, 1): 1.}, 5)
        self.assertRaises(TypeError, eh.set_hopping_def, None, 0, 5)
        self.assertRaises(TypeError, eh.set_hopping_def, None, {0: 1.}, 5)
        self.assertRaises(TypeError, eh.set_hopping_def, None, {(0,): 1.}, 5)
        self.assertRaises(ValueError, eh.set_hopping_def, None, {(0., 1): 1.}, 5)
        self.assertRaises(ValueError, eh.set_hopping_def, None, {(-1, 1): 1.}, 5)
        self.assertRaises(ValueError, eh.set_hopping_def, None, {(0, 0): 1.}, 5)
        self.assertRaises(TypeError, eh.set_hopping_def, None, {(0, 1): 'a'}, 5)

    def test_set_onsite_def(self):
        eh.set_onsite_def({0: 1.}, 5)
        self.assertRaises(TypeError, eh.set_onsite_def, 0, 5)
        self.assertRaises(TypeError, eh.set_onsite_def, {'a': 1.}, 5)
        self.assertRaises(TypeError, eh.set_onsite_def, {0: 'a'}, 5)
        self.assertRaises(ValueError, eh.set_onsite_def, {5: 1.}, 5)

    def test_hop_n1(self):
        hop = np.array([(1,)], dtype=[('n', 'u2')])
        eh.hop_n1(hop)
        empty_hop = np.array([], dtype=[('n', 'u2')])
        self.assertRaises(ValueError, eh.hop_n1, empty_hop)

    def test_empty_onsite(self):
        eh.empty_onsite(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_onsite, np.array([]))

    def test_empty_hop(self):
        eh.empty_hop(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_hop, np.array([]))

    def test_hop_sites(self):
        hop = np.array([(0,), (3,)], dtype=[('i', 'u4')])
        eh.hop_sites(hop, 5)
        self.assertRaises(ValueError, eh.hop_sites, hop, 2)

    def test_empty_ham(self):
        ham_full = sparse.csr_matrix(np.array([[1., 0.], [0., 1.]]))
        eh.empty_ham(ham_full)
        ham_empty = sparse.csr_matrix((2, 2))
        self.assertRaises(RuntimeError, eh.empty_ham, ham_empty)

    def test_empty_en(self):
        eh.empty_en(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_en, np.array([]))

    def test_empty_pola(self):
        eh.empty_pola(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_pola, np.array([]))

    def test_empty_vn(self):
        eh.empty_vn(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_vn, np.array([]))

    def test_empty_ipr(self):
        eh.empty_ipr(np.array([1.]))
        self.assertRaises(RuntimeError, eh.empty_ipr, np.array([]))

    def test_empty_ndarray(self):
        eh.empty_ndarray(np.array([1.]), 'method')
        self.assertRaises(RuntimeError, eh.empty_ndarray, np.array([]), 'method')

    def test_tag(self):
        eh.tag('a', ['a', 'b'])
        self.assertRaises(TypeError, eh.tag, 0, ['a'])
        self.assertRaises(ValueError, eh.tag, 'z', ['a'])

    def test_angle(self):
        eh.angle(45., np.array([45., 135.]), True)
        eh.angle(-135., np.array([45., 135.]), False)
        self.assertRaises(TypeError, eh.angle, 'a', np.array([45.]), True)
        self.assertRaises(ValueError, eh.angle, -45., np.array([45.]), True)
        self.assertRaises(ValueError, eh.angle, 45., np.array([45.]), False)
        self.assertRaises(ValueError, eh.angle, 10., np.array([45.]), True)

    def test_lims(self):
        eh.lims(None)
        eh.lims([1., 2.])
        self.assertRaises(TypeError, eh.lims, 0)
        self.assertRaises(ValueError, eh.lims, (1, 2, 3))
        self.assertRaises(TypeError, eh.lims, ['a', 2.])
        self.assertRaises(ValueError, eh.lims, [2., 1.])

    def test_lims_positive(self):
        eh.lims_positive(None)
        eh.lims_positive([1., 2.])
        self.assertRaises(TypeError, eh.lims_positive, 0)
        self.assertRaises(ValueError, eh.lims_positive, (1, 2, 3))
        self.assertRaises(TypeError, eh.lims_positive, ['a', 2.])
        self.assertRaises(ValueError, eh.lims_positive, [-1., 2.])
        self.assertRaises(ValueError, eh.lims_positive, [2., 1.])


class TestPlotChecks(unittest.TestCase):
    def test_fig(self):
        import matplotlib.pyplot as plt
        fig = plt.figure()
        eh.fig(fig)
        plt.close(fig)
        self.assertRaises(TypeError, eh.fig, 0)

    def test_ani(self):
        class FuncAnimation:
            pass
        eh.ani(FuncAnimation())
        self.assertRaises(TypeError, eh.ani, 0)

    def test_file_format(self):
        eh.file_format('png')
        self.assertRaises(TypeError, eh.file_format, 0)
        self.assertRaises(ValueError, eh.file_format, 'bmp')


class TestPropagationChecks(unittest.TestCase):
    def test_get_pump(self):
        ham_full = sparse.csr_matrix(np.array([[1., 0.], [0., 1.]]))
        eh.get_pump([ham_full])
        self.assertRaises(TypeError, eh.get_pump, 0)
        ham_empty = sparse.csr_matrix((2, 2))
        self.assertRaises(RuntimeError, eh.get_pump, [ham_empty])

    def test_prop_type(self):
        eh.prop_type('real')
        self.assertRaises(TypeError, eh.prop_type, 0)
        self.assertRaises(ValueError, eh.prop_type, 'z')


class TestKSpaceChecks(unittest.TestCase):
    def test_k_vector(self):
        eh.k_vector((1., 2.), 'k', 2)
        eh.k_vector(np.array([1., 2.]), 'k', 2)
        self.assertRaises(TypeError, eh.k_vector, 0, 'k', 2)
        self.assertRaises(ValueError, eh.k_vector, (1.,), 'k', 2)
        self.assertRaises(TypeError, eh.k_vector, (1., 'a'), 'k', 2)

    def test_set_hopping_kspace(self):
        eh.set_hopping_kspace([{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.}], 2, 2)
        self.assertRaises(TypeError, eh.set_hopping_kspace, 0, 2, 2)
        self.assertRaises(TypeError, eh.set_hopping_kspace, [0], 2, 2)
        self.assertRaises(KeyError, eh.set_hopping_kspace, [{'i': 0}], 2, 2)
        self.assertRaises(TypeError, eh.set_hopping_kspace,
                                    [{'i': 0., 'j': 1, 'R': (0, 0), 't': 1.}], 2, 2)
        self.assertRaises(ValueError, eh.set_hopping_kspace,
                                    [{'i': 0, 'j': 5, 'R': (0, 0), 't': 1.}], 2, 2)
        self.assertRaises(ValueError, eh.set_hopping_kspace,
                                    [{'i': 0, 'j': 1, 'R': (0,), 't': 1.}], 2, 2)
        self.assertRaises(TypeError, eh.set_hopping_kspace,
                                    [{'i': 0, 'j': 1, 'R': (0., 0), 't': 1.}], 2, 2)
        self.assertRaises(ValueError, eh.set_hopping_kspace,
                                    [{'i': 0, 'j': 0, 'R': (0, 0), 't': 1.}], 2, 2)
        self.assertRaises(TypeError, eh.set_hopping_kspace,
                                    [{'i': 0, 'j': 1, 'R': (0, 0), 't': 'a'}], 2, 2)

    def test_set_onsite_kspace(self):
        eh.set_onsite_kspace({'a': 1.}, ['a'])
        self.assertRaises(TypeError, eh.set_onsite_kspace, 0, ['a'])
        self.assertRaises(ValueError, eh.set_onsite_kspace, {'z': 1.}, ['a'])
        self.assertRaises(TypeError, eh.set_onsite_kspace, {'a': 'z'}, ['a'])

    def test_k_path_points(self):
        eh.k_path_points([(0., 0.), (1., 1.)], 2)
        self.assertRaises(TypeError, eh.k_path_points, 0, 2)
        self.assertRaises(ValueError, eh.k_path_points, [(0., 0.)], 2)

    def test_nk(self):
        eh.nk(10, 2)
        eh.nk((3, 4), 2)
        self.assertRaises(TypeError, eh.nk, 1.5, 2)
        self.assertRaises(ValueError, eh.nk, 0, 2)
        self.assertRaises(TypeError, eh.nk, [3, 4], 2)
        self.assertRaises(ValueError, eh.nk, (3,), 2)
        self.assertRaises(ValueError, eh.nk, (3, 0), 2)


class TestDosChecks(unittest.TestCase):
    def test_dos_kernel(self):
        eh.dos_kernel('gaussian')
        eh.dos_kernel('lorentzian')
        self.assertRaises(TypeError, eh.dos_kernel, 0)
        self.assertRaises(ValueError, eh.dos_kernel, 'z')


if __name__ == '__main__':
    unittest.main()
