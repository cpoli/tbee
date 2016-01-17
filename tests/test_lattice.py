from TB.lattice import lattice
import unittest
import numpy as np

class TestLattice(unittest.TestCase):
    '''
    Unittest of class **lattice**.
    '''
    def test_unit_cell(self):
        prim_vec = [(0, 1)]
        self.assertRaises(TypeError, lattice, unit_cell=0, prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'z': b'a', 'r0': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'tag': b'a', 'z': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'ab', 'r0': (0, 0)}], prim_vec=prim_vec)
        self.assertRaises(TypeError, lattice, unit_cell=[{'tag': b'a', 'r0': [0, 0]}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'a', 'r0': (0, 'z')}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'a', 'r0': (0, 0, 0)}], prim_vec=prim_vec)

    def test_prim_vec(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec=(0, 1))
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec=[(0)])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0, 0)])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0), (0, 'a')])
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec=[(0, 0), (0, 0), (0, 0)])

    def test_get_lattice(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.get_lattice, n1='z', n2=1)
        self.assertRaises(TypeError, lat.get_lattice, n1=1, n2='z')
        self.assertRaises(ValueError, lat.get_lattice, n1=-1, n2=1)
        self.assertRaises(ValueError, lat.get_lattice, n1=1, n2=-1)
        self.assertRaises(ValueError, lat.get_lattice, n1=10, n2=10)

    def test_add_sites(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.add_sites,
                                    np.array([(0, 0)], dtype=[('x', 'f8'), ('y', 'f8')]))

    def test_remove_sites(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.remove_sites, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.remove_sites, 0)
        self.assertRaises(ValueError, lat.remove_sites, [-1])
        self.assertRaises(ValueError, lat.remove_sites, [10])

    def test_shift_x(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_x, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.shift_x, 0j)

    def test_shift_y(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_y, 5)
        lat.get_lattice(n1=10, n2=1)
        self.assertRaises(TypeError, lat.shift_y, 0j)

    def test_boundary_line(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.boundary_line, cx=0, cy=0, co=0)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.boundary_line, cx=0j, cy=0, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0j, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0, co=0j)

    def test_ellipse(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
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
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__add__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=other)

    def test_iadd(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__iadd__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=other)

    def test_sub(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__sub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=other)

    def test_isub(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__isub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=other)

    def test_plot(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.plot)

    def test_square(self):
        unit_cell = [{'tag': b'a', 'r0': (0, 0)}]
        prim_vec = [(1, 0), (0, 1)]
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1, n2 = 2, 2
        lat.get_lattice(n1=n1, n2=n2)
        lat.shift_x(1.)
        lat.shift_y(1.)
        coor = np.array([(1.0, 1.0, b'a'), (2.0, 1.0, b'a'), (1., 2.0, b'a'), (2.0, 2.0, b'a')], 
                                  dtype=[('x', 'f8'), ('y', 'f8'), ('tag', 'S1')])
        tags = np.array([b'a'])
        sites = 4
        self.assertTrue(np.allclose(lat.coor['x'], coor['x']) == True)
        self.assertTrue(np.allclose(lat.coor['y'], coor['y']) == True)
        self.assertTrue(np.array_equal(lat.coor['tag'], coor['tag']) == True)
        self.assertTrue(np.array_equal(lat.tags, tags) == True)
        self.assertTrue(lat.sites == sites)


if __name__ == '__main__':
    unittest.main()
