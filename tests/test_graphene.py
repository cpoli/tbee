from tbee.graphene import GrapheneLattice, GrapheneSystem, grapheneLat, grapheneSys
import unittest


class TestGrapheneLattice(unittest.TestCase):
    '''
    Unittest of class **GrapheneLattice**.
    '''
    def test_triangle_zigzag(self):
        lat = GrapheneLattice()
        lat.triangle_zigzag(3)
        self.assertTrue(lat.sites > 0)

    def test_hexagon_zigzag(self):
        lat = GrapheneLattice()
        lat.hexagon_zigzag(3)
        self.assertTrue(lat.sites > 0)

    def test_triangle_armchair(self):
        lat = GrapheneLattice()
        lat.triangle_armchair(3)
        self.assertTrue(lat.sites > 0)

    def test_hexagon_armchair(self):
        lat = GrapheneLattice()
        lat.hexagon_armchair(3)
        self.assertTrue(lat.sites > 0)

    def test_square(self):
        lat = GrapheneLattice()
        lat.square(3)
        self.assertTrue(lat.sites > 0)

    def test_circle_even(self):
        lat = GrapheneLattice()
        lat.circle(4)
        self.assertTrue(lat.sites > 0)

    def test_circle_odd(self):
        lat = GrapheneLattice()
        lat.circle(3)
        self.assertTrue(lat.sites > 0)


class TestGrapheneSystem(unittest.TestCase):
    '''
    Unittest of class **GrapheneSystem**.
    '''
    def test_strain_and_butterfly(self):
        lat = GrapheneLattice()
        lat.hexagon_zigzag(3)
        sys = GrapheneSystem(lat)
        sys.set_hop_linear_strain(t=1., beta=0.)
        sys.get_ham()
        sys.get_eig()
        self.assertTrue(len(sys.en) == lat.sites)
        beta_lims = sys.get_beta_lims()
        self.assertEqual(len(beta_lims), 2)
        sys.get_butterfly(t=1., N=3)
        self.assertEqual(sys.butterfly.shape, (3, lat.sites))
        self.assertEqual(len(sys.betas), 3)

    def test_backward_compatible_aliases(self):
        self.assertIs(grapheneLat, GrapheneLattice)
        self.assertIs(grapheneSys, GrapheneSystem)


if __name__ == '__main__':
    unittest.main()
