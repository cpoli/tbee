from tbee.lattice import lattice
import unittest
import numpy as np

class TestLattice(unittest.TestCase):
    '''
    Unittest of class **lattice**.
    '''
    def test_unit_cell(self):
        prim_vec = [(0, 1)]
        self.assertRaises(TypeError, lattice, unit_cell=0, prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'z': 'a', 'r0': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'tag': 'a', 'z': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': 'ab', 'r0': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(TypeError, lattice, unit_cell=[{'tag': 'a', 'r0': [0, 0]}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': 'a', 'r0': (0, 'z')}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': 'a', 'r0': (0, 0, 0)}], prim_vec=prim_vec)

    def test_prim_vec(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec=(0, 1))
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec=[(0)])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0, 0)])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0), (0, 'a')])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0), (0, 0), (0, 0)])

    def test_get_lattice(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.get_lattice, n1='z', n2=1)
        self.assertRaises(TypeError, lat.get_lattice, n1=1, n2='z')
        self.assertRaises(ValueError, lat.get_lattice, n1=-1, n2=1)
        self.assertRaises(ValueError, lat.get_lattice, n1=1, n2=-1)
        self.assertRaises(ValueError, lat.get_lattice, n1=10, n2=10)

    def test_add_sites(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.add_sites,
                                    np.array([(0, 0)], dtype=[('x', 'f8'), ('y', 'f8')]))

    def test_remove_sites(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.remove_sites, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.remove_sites, 0)
        self.assertRaises(ValueError, lat.remove_sites, [-1])
        self.assertRaises(ValueError, lat.remove_sites, [10])

    def test_shift_x(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_x, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.shift_x, 0j)

    def test_shift_y(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_y, 5)
        lat.get_lattice(n1=10, n2=1)
        self.assertRaises(TypeError, lat.shift_y, 0j)

    def test_boundary_line(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.boundary_line, cx=0, cy=0, co=0)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.boundary_line, cx=0j, cy=0, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0j, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0, co=0j)

    def test_ellipse(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.ellipse_in, rx=2, ry=2, x0=0., y0=0.)
        self.assertRaises(RuntimeError, lat.ellipse_out, rx=2, ry=2, x0=0., y0=0.)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.ellipse_in, rx=2j, ry=2, x0=0., y0=0.)
        self.assertRaises(TypeError, lat.ellipse_in, rx=2, ry=2j, x0=0., y0=0.)
        self.assertRaises(ValueError, lat.ellipse_in, rx=0, ry=2, x0=0., y0=0.)
        self.assertRaises(ValueError, lat.ellipse_in, rx=2, ry=0, x0=0., y0=0.)
        self.assertRaises(TypeError, lat.ellipse_out, rx=2j, ry=2, x0=0., y0=0.)
        self.assertRaises(TypeError, lat.ellipse_out, rx=2, ry=2j, x0=0., y0=0.)
        self.assertRaises(ValueError, lat.ellipse_out, rx=0, ry=2, x0=0., y0=0.)
        self.assertRaises(ValueError, lat.ellipse_out, rx=2, ry=0, x0=0., y0=0.)

    def test_add(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__add__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=other)

    def test_iadd(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__iadd__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=other)

    def test_sub(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__sub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=other)

    def test_isub(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__isub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=other)

    def test_plot(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.plot)

    def test_square(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1, n2 = 2, 2
        lat.get_lattice(n1=n1, n2=n2)
        lat.shift_x(1.)
        lat.shift_y(1.)
        coor = np.array([(1.0, 1.0, 'a'), (2.0, 1.0, 'a'), (1., 2.0, 'a'), (2.0, 2.0, 'a')],
                                  dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'U1')])
        tags = np.array(['a'])
        sites = 4
        self.assertTrue(np.allclose(lat.coor['x'], coor['x']) == True)
        self.assertTrue(np.allclose(lat.coor['y'], coor['y']) == True)
        self.assertTrue(np.array_equal(lat.coor['tag'], coor['tag']) == True)
        self.assertTrue(np.array_equal(lat.tags, tags) == True)
        self.assertTrue(lat.sites == sites)

    def test_add_sites_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat.get_lattice(n1=2, n2=2)
        sites_before = lat.sites
        coor = np.array([(-1., -1., 'b'), (-2., -2., 'c')],
                                  dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'U1')])
        lat.add_sites(coor)
        self.assertEqual(lat.sites, sites_before + 2)
        self.assertTrue(set('b c a'.split()) <= set(lat.tags.tolist()))

    def test_remove_dangling(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.remove_dangling)
        lat.get_lattice(n1=5, n2=5)
        lat.remove_dangling()
        self.assertEqual(lat.sites, len(lat.coor))

    def test_change_sign(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.change_sign_x)
        self.assertRaises(RuntimeError, lat.change_sign_y)
        lat.get_lattice(n1=3, n2=3)
        x, y = lat.coor['x'].copy(), lat.coor['y'].copy()
        lat.change_sign_x()
        lat.change_sign_y()
        self.assertTrue(np.allclose(lat.coor['x'], -x))
        self.assertTrue(np.allclose(lat.coor['y'], -y))

    def test_ellipse_out_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat.get_lattice(n1=10, n2=10)
        sites_before = lat.sites
        lat.ellipse_out(rx=1., ry=1., x0=4.5, y0=4.5)
        self.assertTrue(lat.sites < sites_before)

    def test_center(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.center)
        lat.get_lattice(n1=4, n2=4)
        lat.center()
        self.assertTrue(np.isclose(np.mean(lat.coor['x']), 0.))
        self.assertTrue(np.isclose(np.mean(lat.coor['y']), 0.))

    def test_rotation(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.rotation, 90.)
        lat.get_lattice(n1=4, n2=4)
        self.assertRaises(TypeError, lat.rotation, '90')
        lat.rotation(360.)
        self.assertTrue(lat.sites == 16)

    def test_clean_coor(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.clean_coor)
        lat.get_lattice(n1=3, n2=3)
        duplicate = lat.coor[:1].copy()
        lat.add_sites(duplicate)
        self.assertEqual(lat.sites, 10)
        lat.clean_coor()
        self.assertEqual(lat.sites, 9)

    def test_add_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat1.get_lattice(n1=2, n2=2)
        lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat2.get_lattice(n1=2, n2=2)
        lat2.shift_x(10.)
        lat3 = lat1 + lat2
        self.assertEqual(lat3.sites, lat1.sites + lat2.sites)

    def test_iadd_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat1.get_lattice(n1=2, n2=2)
        sites_before = lat1.sites
        lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat2.get_lattice(n1=2, n2=2)
        lat2.shift_x(10.)
        lat1 += lat2
        self.assertEqual(lat1.sites, sites_before + lat2.sites)

    def test_sub_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat1.get_lattice(n1=3, n2=3)
        lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat2.get_lattice(n1=2, n2=2)
        lat3 = lat1 - lat2
        self.assertEqual(lat3.sites, lat1.sites - lat2.sites)

    def test_isub_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat1.get_lattice(n1=3, n2=3)
        sites_before = lat1.sites
        lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat2.get_lattice(n1=2, n2=2)
        lat1 -= lat2
        self.assertEqual(lat1.sites, sites_before - lat2.sites)

    def test_plot_valid(self):
        unit_cell = [{'tag': 'a', 'r0': (0, 0)}, {'tag': 'b', 'r0': (0.5, 0.5)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat.get_lattice(n1=3, n2=3)
        fig = lat.plot(plt_index=True, axis=True, figsize=(4, 4))
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = lat.plot()
        self.assertEqual(fig2.__class__.__name__, 'Figure')
        lat.show()


if __name__ == '__main__':
    unittest.main()
