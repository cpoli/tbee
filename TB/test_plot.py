from lattice import *
from system import *
from plot import *

import unittest
import numpy as np


def init():
    unit_cell = [{'tag': b'a', 'r0': [0, 0]}]
    prim_vec = {'norm': 1, 'angle': 0}
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=5, n2=5)
    sys = system(lat)
    sys.set_hopping([{'n': 1,'t':1.}])
    sys.get_ham()
    sys.get_eig(eigenvec=True)
    return lat, sys


class TestPlot(unittest.TestCase):
    '''
    Unittest of class **plot**.
    '''
    def test_plot(self):
        self.assertRaises(TypeError, plot, lat=None, sys=None)
        self.assertRaises(TypeError, plot, lat=None, sys=0)
        self.assertRaises(TypeError, plot, lat=0, sys=None)

    def test_lattice(self):
        lat, sys = init()
        plt = plot(lat=lat)
        self.assertRaises(TypeError, plt.lattice, ms=0j)
        self.assertRaises(ValueError, plt.lattice, ms=0)
        self.assertRaises(TypeError, plt.lattice, lw=0j)
        self.assertRaises(ValueError, plt.lattice, lw=0)
        self.assertRaises(TypeError, plt.lattice, fs=0j)
        self.assertRaises(ValueError, plt.lattice, fs=0)
        self.assertRaises(TypeError, plt.lattice, plt_hop='a')
        self.assertRaises(TypeError, plt.lattice, plt_index='a')
        self.assertRaises(TypeError, plt.lattice, tics='a')
        self.assertRaises(TypeError, plt.lattice, figsize=0)
        self.assertRaises(ValueError, plt.lattice, figsize=[0, 0, 0])
        self.assertRaises(ValueError, plt.lattice, figsize=[2, -1])
        self.assertRaises(ValueError, plt.lattice, figsize=[-1, 2])


if __name__ == '__main__':
    unittest.main()