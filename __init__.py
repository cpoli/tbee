# Copyright 2014 Charles Poli.
#
# This file is part of TB.  It is subject to the license terms in the
# LICENSE file found in the top-level directory of this distribution and at
# https://github.com/cpoli/TB.  

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


__all__ = ['latticeTB', 'eigTB', 'plotTB', 'propTB', 'grapheneTB', 'liebTB']
from latticeTB import latticeTB
from eigTB import eigTB
from plotTB import plotTB
from propTB import propTB
from grapheneTB import grapheneTB
from liebTB import liebTB

#    exec('from . import {}'.format(module))

#execfile(test_latticeTB.py)
#execfile(test_eigTB.py, "-v")
