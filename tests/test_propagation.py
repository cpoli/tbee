from lattice import *
from system import *
from propagation import *
import unittest
import numpy as np
from math import sqrt


PI = np.pi


def init():
    unit_cell = [{'tag': b'a', 'r0': (0., 0.)}, 
                      {'tag': b'b', 'r0': (0.5, 0.5/sqrt(3))}]
    prim_vec = [(1, 0.), 
                (cos(PI/3), sin(PI/3))]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    n1, n2 = 2, 2
    lat.get_lattice(n1=n1, n2=n2)
    lat.remove_dangling()
    sys = system(lat=lat)
    sys.set_hopping([{'n': 1, 't': 1.}])
    sys.get_ham()
    return sys