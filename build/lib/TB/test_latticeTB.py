from latticeTB import latticeTB
import unittest
import numpy as np

class TestLatticeTB(unittest.TestCase):
    '''
    Unittest class of **latticeTB**.
    '''

    def test_errors_ri(self):
        # ri not a list.
        self.assertRaises(TypeError, latticeTB, None, 0, 0., 0.)
        # ri is not a container of 2D coordinates.
        self.assertRaises(ValueError, latticeTB, [], [], 0., 0.)
        self.assertRaises(ValueError, latticeTB, [[0.]], 0., 0., 0.) 
        self.assertRaises(ValueError, latticeTB, [[0., 0., 0.]], 0, 0., 0.) 

    def test_errors_tags(self):
        # tags not a container of binary chars.
        self.assertRaises(TypeError, latticeTB, [[0., 0.]], None, 1., 1.)
        self.assertRaises(ValueError, latticeTB, [[0., 0.]], ['a'], 1., 1.)
        # ri and tags have Different lengths 
        self.assertRaises(ValueError, latticeTB, [[0., 0.], [1., 1.]], [b'a'], 1., 1.)

    def test_errors_ri_tags(self):
        # ri and tags not of same length.
        self.assertRaises(ValueError, latticeTB, [[0., 0.], [1., 1.]], [b'a'], 1., 1.)

    def test_errors_nor(self):
        # nor not a positive number
        self.assertRaises(TypeError, latticeTB, [[0., 0.]],  [b'a'], 1+1j, 1.)

    def test_errors_ang(self):
        # nor not a real number
        self.assertRaises(TypeError, latticeTB, [[0., 0.]],  [b'a'], 1., 1+1j)

    def test_remove_sites(self):
        tb = latticeTB([[0., 0.]],  [b'a'], 1.,0.)
        # method get_lattice not called
        self.assertRaises(RuntimeError, tb.remove_sites, 5)
        tb.get_lattice(nx=10, ny=1)
        # input parameter not a list
        self.assertRaises(TypeError, tb.remove_sites, 0)
        # input parameter contains negative numbers
        self.assertRaises(IndexError, tb.remove_sites, [-1])
        # input parameter contains integers larger than sites-1
        self.assertRaises(IndexError, tb.remove_sites, [10])

    def test_remove_dangling(self):
        tb = latticeTB([[0., 0.], [1., 0.], [0., 1.]],  [b'a', b'b', b'c'], 2., 0.)
        # method get_lattice not called
        self.assertRaises(RuntimeError, tb.remove_dangling, 1)
        tb.get_lattice(nx=2, ny=2)
        # method get_lattice not called
        self.assertRaises(TypeError, tb.remove_dangling, None)
        self.assertRaises(TypeError, tb.remove_dangling, -1)

    def test_plt_lattice(self):
        tb = latticeTB([[0., 0.], [1., 0.], [0., 1.]],  [b'a', b'b', b'c'], 2., 0.)
        self.assertRaises(RuntimeError, tb.plt_lattice)
        tb.get_lattice(nx=2, ny=2)
        self.assertRaises(TypeError, tb.plt_lattice, ms=-1)
        self.assertRaises(TypeError, tb.plt_lattice, fs=-1)
        self.assertRaises(TypeError, tb.plt_lattice, plt_index=None)
        self.assertRaises(ValueError, tb.plt_lattice, figsize=[-1, 1])
        self.assertRaises(ValueError, tb.plt_lattice, figsize=[0, 1, 2])
        self.assertRaises(ValueError, tb.plt_lattice, colors=['r'])
        self.assertRaises(ValueError, tb.plt_lattice, colors=['a', 'b', 'c'])
        self.assertRaises(ValueError, tb.plt_lattice, colors=['#008000', '#00FFFF', '#0000F'])


if __name__ == '__main__':
    unittest.main()
