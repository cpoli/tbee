from latticeTB import latticeTB
import unittest
import numpy as np

class TestLatticeTB(unittest.TestCase):
    '''
    Unittest class of **latticeTB**.
    '''
    def test_unit_cell(self):
        prim_vec = {'norm': 1, 'angle': 0}
        self.assertRaises(TypeError, latticeTB, unit_cell=0, prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(KeyError, latticeTB, unit_cell=[{'z': b'a', 'r0': [0, 0]}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(KeyError, latticeTB, unit_cell=[{'tag': b'a', 'z': [0, 0]}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(TypeError, latticeTB, unit_cell=[{'tag': 'ab', 'r0': [0, 0]}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(ValueError, latticeTB, unit_cell=[{'tag': b'ab', 'r0': [0, 0]}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(TypeError, latticeTB, unit_cell=[{'tag': b'a', 'r0': 0}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(ValueError, latticeTB, unit_cell=[{'tag': b'a', 'r0': [0, 'z']}], prim_vec={'angle':0, 'norm': 1})
        self.assertRaises(ValueError, latticeTB, unit_cell=[{'tag': b'a', 'r0': [0, 0, 0]}], prim_vec={'angle':0, 'norm': 1})

    def test_prim_vec(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        self.assertRaises(TypeError, latticeTB, unit_cell=unit_cell, prim_vec=[])
        self.assertRaises(KeyError, latticeTB, unit_cell=unit_cell, prim_vec={'z':0, 'norm': 1})
        self.assertRaises(KeyError, latticeTB, unit_cell=unit_cell, prim_vec={'angle':0, 'z': 1})
        self.assertRaises(TypeError, latticeTB, unit_cell=unit_cell, prim_vec={'angle':'z', 'norm': 1})
        self.assertRaises(TypeError, latticeTB, unit_cell=unit_cell, prim_vec={'angle':0, 'norm': 'z'})
        self.assertRaises(TypeError, latticeTB, unit_cell=unit_cell, prim_vec={'angle':0j, 'norm': -1})
        self.assertRaises(ValueError, latticeTB, unit_cell=unit_cell, prim_vec={'angle':0, 'norm': -1})

    def test_get_lattice(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(TypeError, lat.get_lattice, n1='z', n2=1)
        self.assertRaises(TypeError, lat.get_lattice, n1=1, n2='z')
        self.assertRaises(ValueError, lat.get_lattice, n1=-1, n2=1)
        self.assertRaises(ValueError, lat.get_lattice, n1=1, n2=-1)

    def test_remove_sites(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.remove_sites, 5)
        lat.get_lattice(n1=10)
        self.assertRaises(TypeError, lat.remove_sites, 0)
        self.assertRaises(ValueError, lat.remove_sites, [-1])
        self.assertRaises(ValueError, lat.remove_sites, [10])
        
    def test_coor_cut_x(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 90}
        lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.coor_cut_x, 5)
        lat.get_lattice(n1=10, n2=10)
        self.assertRaises(TypeError, lat.coor_cut_x, 0)
        self.assertRaises(TypeError, lat.coor_cut_x, ['z', 5])
        self.assertRaises(TypeError, lat.coor_cut_x, [5, 'z'])
        self.assertRaises(ValueError, lat.coor_cut_x, [5, 0])
    '''
    def test_plt_lattice(self):
        unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
        prim_vec = {'norm': 1, 'angle': 0}
        lat = latticeTB(unit_cell=unit_cell, prim_vec=prim_vec)
        self.assertRaises(RuntimeError, lat.plt_lattice)
        lat.get_lattice(nx=2, ny=2)
        self.assertRaises(TypeError, lat.plt_lattice, ms=-1)
        self.assertRaises(TypeError, lat.plt_lattice, fs=-1)
        self.assertRaises(TypeError, lat.plt_lattice, plt_index=None)
        self.assertRaises(ValueError, lat.plt_lattice, figsize=[-1, 1])
        self.assertRaises(ValueError, lat.plt_lattice, figsize=[0, 1, 2])
        self.assertRaises(ValueError, lat.plt_lattice, colors=['r'])
        self.assertRaises(ValueError, lat.plt_lattice, colors=['a', 'b', 'c'])
        self.assertRaises(ValueError, lat.plt_lattice, colors=['#008000', '#00FFFF', '#0000F'])
    '''

if __name__ == '__main__':
    unittest.main()
