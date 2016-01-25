# Copyright 2014 Charles Poli.
#
# This file is part of TBEE.  It is subject to the license terms in the
# LICENSE file found in the top-level directory of this distribution and at
# https://github.com/cpoli/tbee.  

import numpy as np
import numpy.core.defchararray as npc
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.legend_handler import HandlerLine2D
from math import sqrt, pi, sin, cos
PI = pi


from distutils.core import setup

setup(name='tbee',
      version='0.1',
      py_modules=['lattice', 'system', 'plot', 'propagation', 'save', 'error_handling'],
      )

__all__ = ["lattice", "system", "plot", "propagation", "save", "error_handling"]

from tbee.lattice import *
from tbee.system import *
from tbee.plot import *
from tbee.propagation import *
from tbee.save import *
import tbee.error_handling
