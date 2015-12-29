from lattice import lattice
import unittest
import numpy as np

class TestLattice(unittest.TestCase):
    '''
    Unittest of class **lattice**.
    '''
    def test_unit_cell(self):
        prim_vec = {'angle':0, 'norm': 1}
        self.assertRaises(TypeError, lattice, unit_cell=0, prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'z': b'a', 'r0': [0, 0]}], prim_vec=prim_vec)
        self.assertRaises(KeyError, lattice, unit_cell=[{'tag': b'a', 'z': [0, 0]}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'ab', 'r0': [0, 0]}], prim_vec=prim_vec)
        self.assertRaises(TypeError, lattice, unit_cell=[{'tag': b'a', 'r0': 0}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'a', 'r0': [0, 'z']}], prim_vec=prim_vec)
        self.assertRaises(ValueError, lattice, unit_cell=[{'tag': b'a', 'r0': [0, 0, 0]}], prim_vec=prim_vec)

    def test_prim_vec(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec=[])
        self.assertRaises(KeyError, lattice, unit_cell=unit_cell, prim_vec={'z':0, 'norm': 1})
        self.assertRaises(KeyError, lattice, unit_cell=unit_cell, prim_vec={'angle':0, 'z': 1})
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec={'angle':'z', 'norm': 1})
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec={'angle':0, 'norm': 'z'})
        self.assertRaises(TypeError, lattice, unit_cell=unit_cell, prim_vec={'angle':0j, 'norm': -1})
        self.assertRaises(ValueError, lattice, unit_cell=unit_cell, prim_vec={'angle':0, 'norm': -1})

    def test_get_lattice(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.get_lattice, n1='z', n2=1)
        self.assertRaises(TypeError, lat.get_lattice, n1=1, n2='z')
        self.assertRaises(ValueError, lat.get_lattice, n1=-1, n2=1)
        self.assertRaises(ValueError, lat.get_lattice, n1=1, n2=-1)

    def test_add_sites(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.add_sites,
                                    np.array([(0, 0)], dtype=[('x', 'f16'), ('y', 'f16')]))

    def test_remove_sites(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.remove_sites, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.remove_sites, 0)
        self.assertRaises(ValueError, lat.remove_sites, [-1])
        self.assertRaises(ValueError, lat.remove_sites, [10])
        
    def test_shift_x(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_x, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.shift_x, 0j)

    def test_shift_y(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.shift_x, 5)
        lat.get_lattice(n1=1, n2=10)
        self.assertRaises(TypeError, lat.shift_x, 0j)

    def test_boundary_line(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.boundary_line, cx=0, cy=0, co=0)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.boundary_line, cx=0j, cy=0, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0j, co=0)
        self.assertRaises(TypeError, lat.boundary_line, cx=0, cy=0, co=0j)

    def test_ellipse(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.ellipse_in, a=2, b=2)
        self.assertRaises(RuntimeError, lat.ellipse_out, a=2, b=2)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.ellipse_in, a=2j, b=2)
        self.assertRaises(TypeError, lat.ellipse_in, a=2, b=2j)
        self.assertRaises(ValueError, lat.ellipse_in, a=0, b=2)
        self.assertRaises(ValueError, lat.ellipse_in, a=2, b=0)
        self.assertRaises(TypeError, lat.ellipse_out, a=2j, b=2)
        self.assertRaises(TypeError, lat.ellipse_out, a=2, b=2j)
        self.assertRaises(ValueError, lat.ellipse_out, a=0, b=2)
        self.assertRaises(ValueError, lat.ellipse_out, a=2, b=0)

    def test_add(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__add__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__add__, other=other)

    def test_iadd(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__iadd__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__iadd__, other=other)

    def test_isub(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__sub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__sub__, other=other)

    def test_isub(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=lat)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.__isub__, other=0)
        other = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.__isub__, other=other)

    def test_square(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1, n2 = 2, 2
        lat.get_lattice(n1=n1, n2=n2)
        lat.shift_x(1.)
        lat.shift_y(1.)
        coor = np.array([(1.0, 1.0, b'a'), (2.0, 1.0, b'a'), (1., 2.0, b'a'), (2.0, 2.0, b'a')], 
                                  dtype=[('x', 'f16'), ('y', 'f16'), ('tag', 'S1')])
        tags = np.array([b'a'])
        sites = 4
        self.assertTrue(np.allclose(lat.coor['x'], coor['x']) == True)
        self.assertTrue(np.allclose(lat.coor['y'], coor['y']) == True)
        self.assertTrue(np.array_equal(lat.coor['tag'], coor['tag']) == True)
        self.assertTrue(np.array_equal(lat.tags, tags) == True)
        self.assertTrue(lat.sites == sites)


    def test_graphene(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}, 
                          {'tag': b'b', 'r0': [0.5*np.sqrt(3), 0.5]}]
        prim_vec = {'norm': np.sqrt(3), 'angle': 60}
        lat1 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        n1, n2 = 2, 2
        lat1.get_lattice(n1=n1, n2=n2)
        lat2 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat2.get_lattice(n1=n1, n2=n2)
        lat2.shift_x(2*np.sqrt(3))
        lat2.shift_y(3)
        lat3 = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        lat3.get_lattice(n1=n1, n2=n2)
        lat3.shift_x(4*np.sqrt(3))
        lat3.shift_y(6)
        lat = lat3 + lat2 + lat1
        tags = np.array([b'a', b'b'])
        sites = 8
        coor = np.array([(6.9282032302755087727, 6.0, b'a'),
                                 (8.6602540378443859659, 6.0, b'a'),
                                 (7.7942286340599473693, 6.5, b'b'),
                                 (9.5262794416288245625, 6.5, b'b'),
                                 (7.7942286340599475636, 7.499999999999999778, b'a'),
                                 (9.5262794416288247568, 7.499999999999999778, b'a'),
                                 (8.6602540378443861602, 7.999999999999999778, b'b'),
                                 (10.392304845413263353, 7.999999999999999778, b'b'),
                                 (3.4641016151377543864, 3.0, b'a'),
                                 (5.1961524227066315795, 3.0, b'a'),
                                 (4.3301270189221929829, 3.5, b'b'),
                                 (6.0621778264910701761, 3.5, b'b'),
                                 (4.3301270189221931772, 4.499999999999999778, b'a'),
                                 (6.0621778264910703704, 4.499999999999999778, b'a'),
                                 (5.1961524227066317738, 4.999999999999999778, b'b'),
                                 (6.928203230275508967, 4.999999999999999778, b'b'),
                                 (0.0, 0.0, b'a'), (1.7320508075688771932, 0.0, b'a'),
                                 (0.86602540378443859659, 0.5, b'b'),
                                 (2.5980762113533157898, 0.5, b'b'),
                                 (0.86602540378443879077, 1.499999999999999778, b'a'),
                                 (2.5980762113533159841, 1.499999999999999778, b'a'),
                                 (1.7320508075688773874, 1.999999999999999778, b'b'),
                                 (3.4641016151377545806, 1.999999999999999778, b'b')], 
                                dtype=[('x', '<f16'), ('y', '<f16'), ('tag', 'S1')])
        self.assertTrue(np.allclose(lat.coor['x'], coor['x']) == True)
        self.assertTrue(np.allclose(lat.coor['y'], coor['y']) == True)
        self.assertTrue(np.array_equal(lat.coor['tag'], coor['tag']) == True)
        self.assertTrue(np.array_equal(lat.tags, tags) == True)
        self.assertTrue(lat.sites == sites)


if __name__ == '__main__':
    unittest.main()
