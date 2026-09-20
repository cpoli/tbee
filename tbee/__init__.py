# Copyright 2014 Charles Poli.
#
# This file is part of TBEE.  It is subject to the license terms in the
# LICENSE file found in the top-level directory of this distribution and at
# https://github.com/cpoli/tbee.

"""tbee: build and solve Tight-Binding models."""

__version__ = "0.2.0"

__all__ = [
    "Lattice", "System", "Plot", "Propagation", "Save", "KSpace",
    "reciprocal_vectors", "error_handling",
]

# NOTE: these are explicit imports, not `from tbee.<module> import *`.
# A wildcard import here would rebind the `tbee.<module>` submodule
# attributes to the classes they define (since e.g. tbee/lattice.py both
# *is* the submodule `tbee.lattice` and defines a `lattice` alias of the
# same name), breaking `import tbee.lattice as lattice`-style imports.
from tbee.lattice import Lattice
from tbee.system import System
from tbee.plot import Plot
from tbee.propagation import Propagation
from tbee.save import Save
from tbee.kspace import KSpace, reciprocal_vectors
import tbee.error_handling
