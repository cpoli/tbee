import unittest
import numpy as np

import tbee.lattices as lattices
from tbee.kspace import KSpace, reciprocal_vectors


def bz_grid(prim_vec, n=20):
    b1, b2 = (np.array(v) for v in reciprocal_vectors(prim_vec))
    return np.array([[i/n*b1[0] + j/n*b2[0], i/n*b1[1] + j/n*b2[1]]
                            for i in range(n) for j in range(n)])


class TestLattices(unittest.TestCase):
    '''
    Unittest of module **lattices**.
    '''
    def test_errors(self):
        for fn in (lattices.chain, lattices.square, lattices.triangular,
                          lattices.honeycomb, lattices.kagome, lattices.lieb):
            self.assertRaises(TypeError, fn, 'a')
            self.assertRaises(ValueError, fn, -1.)

    def test_chain(self):
        lat = lattices.chain(a=2.)
        self.assertEqual(len(lat.unit_cell), 1)
        lat.get_lattice(n1=5)
        d = np.diff(lat.coor['x'])
        self.assertTrue(np.allclose(d, 2.))

    def test_square(self):
        lat = lattices.square(a=1.5)
        self.assertEqual(len(lat.unit_cell), 1)
        lat.get_lattice(n1=3, n2=3)
        self.assertEqual(lat.sites, 9)

    def test_triangular_coordination(self):
        lat = lattices.triangular()
        lat.get_lattice(n1=5, n2=5)
        lat.center()
        dif_x = lat.coor['x'] - lat.coor['x'].reshape(-1, 1)
        dif_y = lat.coor['y'] - lat.coor['y'].reshape(-1, 1)
        dist = np.sqrt(dif_x**2 + dif_y**2)
        # a bulk site of a triangular lattice has 6 nearest neighbors.
        i_bulk = np.argmin(dist.max(axis=1))
        n_nn = np.sum(np.isclose(dist[i_bulk], 1., atol=1e-6))
        self.assertEqual(n_nn, 6)

    def test_honeycomb_matches_graphene_convention(self):
        from tbee.graphene import GrapheneLattice
        lat = lattices.honeycomb(a=1.)
        gra = GrapheneLattice()
        self.assertEqual(lat.unit_cell, gra.unit_cell)
        self.assertEqual(lat.prim_vec, gra.prim_vec)

    def test_kagome_flat_band(self):
        lat = lattices.kagome()
        kag = KSpace(lat)
        t = 1.
        kag.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t},
                                {'i': 0, 'j': 1, 'R': (-1, 0), 't': t},
                                {'i': 0, 'j': 2, 'R': (0, 0), 't': t},
                                {'i': 0, 'j': 2, 'R': (0, -1), 't': t},
                                {'i': 1, 'j': 2, 'R': (0, 0), 't': t},
                                {'i': 1, 'j': 2, 'R': (1, -1), 't': t}])
        self.assertEqual(len(kag._hop), 12)  # 4 neighbors x 3 sites
        en = kag.get_bands(bz_grid(lat.prim_vec))
        self.assertTrue(np.allclose(en[:, 0], -2*t, atol=1e-8))
        self.assertFalse(np.allclose(en[:, 1], en[0, 1], atol=1e-6))

    def test_lieb_flat_band(self):
        lat = lattices.lieb()
        lb = KSpace(lat)
        t = 1.
        lb.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t},
                              {'i': 0, 'j': 1, 'R': (-1, 0), 't': t},
                              {'i': 0, 'j': 2, 'R': (0, 0), 't': t},
                              {'i': 0, 'j': 2, 'R': (0, -1), 't': t}])
        en = lb.get_bands(bz_grid(lat.prim_vec))
        self.assertTrue(np.allclose(en[:, 1], 0., atol=1e-8))
        self.assertTrue(np.allclose(np.sort(en[0]), [-2*np.sqrt(2)*t, 0., 2*np.sqrt(2)*t]))


if __name__ == '__main__':
    unittest.main()
