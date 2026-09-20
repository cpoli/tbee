"""
A small library of common 2D Bravais lattices, ready to feed into
:class:`tbee.lattice.Lattice`, :class:`tbee.system.System`, or
:class:`tbee.kspace.KSpace`.

Each function returns a fresh ``Lattice`` instance with *unit_cell* and
*prim_vec* already set (call ``get_lattice`` yourself to build a finite
flake, or hand it straight to ``KSpace`` for a periodic/band-structure
calculation). Nearest-neighbor sites are a distance *a* apart.

Example usage::

    import tbee.lattices as lattices
    lat = lattices.kagome()
    lat.get_lattice(n1=6, n2=6)
"""
from __future__ import annotations

from math import sqrt

from tbee.lattice import Lattice
import tbee.error_handling as error_handling


def _lat(unit_cell: list[dict], prim_vec: list[tuple[float, float]]) -> Lattice:
    return Lattice(unit_cell=unit_cell, prim_vec=prim_vec)


def chain(a: float = 1.) -> Lattice:
    '''
    1D chain: one site per unit cell.

    :param a: Positive real number. Default value 1. Lattice constant.
    '''
    error_handling.positive_real(a, 'a')
    return _lat([{'tag': 'a', 'r0': (0., 0.)}], [(a, 0.)])


def square(a: float = 1.) -> Lattice:
    '''
    Square lattice: one site per unit cell.

    :param a: Positive real number. Default value 1. Lattice constant.
    '''
    error_handling.positive_real(a, 'a')
    return _lat([{'tag': 'a', 'r0': (0., 0.)}], [(a, 0.), (0., a)])


def triangular(a: float = 1.) -> Lattice:
    '''
    Triangular lattice: one site per unit cell.

    :param a: Positive real number. Default value 1. Lattice constant.
    '''
    error_handling.positive_real(a, 'a')
    return _lat([{'tag': 'a', 'r0': (0., 0.)}],
                       [(a, 0.), (0.5*a, 0.5*sqrt(3)*a)])


def honeycomb(a: float = 1.) -> Lattice:
    '''
    Honeycomb lattice (e.g. graphene): two sites per unit cell, nearest
    neighbors a distance *a* apart. See also :class:`tbee.graphene.GrapheneLattice`
    for ready-made finite flakes of various shapes.

    :param a: Positive real number. Default value 1. Nearest-neighbor distance.
    '''
    error_handling.positive_real(a, 'a')
    dx, dy = 0.5*sqrt(3)*a, 0.5*a
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (dx, dy)}]
    prim_vec = [(2*dx, 0.), (dx, 1.5*a)]
    return _lat(unit_cell, prim_vec)


def kagome(a: float = 1.) -> Lattice:
    '''
    Kagome lattice: three sites per unit cell (tags 'a', 'b', 'c'), arranged
    as corner-sharing triangles on a triangular Bravais lattice. With
    uniform nearest-neighbor hopping, this lattice famously has an exactly
    flat band (at E = -2t for hopping amplitude t).

    :param a: Positive real number. Default value 1. Nearest-neighbor distance.

    Example usage (nearest-neighbor hopping, for :class:`tbee.kspace.KSpace`)::

        lat = lattices.kagome()
        kag = KSpace(lat)
        kag.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t},
                                {'i': 0, 'j': 1, 'R': (-1, 0), 't': t},
                                {'i': 0, 'j': 2, 'R': (0, 0), 't': t},
                                {'i': 0, 'j': 2, 'R': (0, -1), 't': t},
                                {'i': 1, 'j': 2, 'R': (0, 0), 't': t},
                                {'i': 1, 'j': 2, 'R': (1, -1), 't': t}])
    '''
    error_handling.positive_real(a, 'a')
    a1 = (2*a, 0.)
    a2 = (a, sqrt(3)*a)
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)},
                       {'tag': 'b', 'r0': (a, 0.)},
                       {'tag': 'c', 'r0': (0.5*a, 0.5*sqrt(3)*a)}]
    return _lat(unit_cell, [a1, a2])


def lieb(a: float = 1.) -> Lattice:
    '''
    Lieb lattice: three sites per unit cell (tag 'a': corner site; tags
    'b', 'c': edge-center sites) on a square Bravais lattice. With uniform
    nearest-neighbor hopping, this lattice famously has an exactly flat
    band (at E = 0), squeezed between two dispersive bands.

    :param a: Positive real number. Default value 1. Nearest-neighbor distance.

    Example usage (nearest-neighbor hopping, for :class:`tbee.kspace.KSpace`)::

        lat = lattices.lieb()
        lb = KSpace(lat)
        lb.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t},
                               {'i': 0, 'j': 1, 'R': (-1, 0), 't': t},
                               {'i': 0, 'j': 2, 'R': (0, 0), 't': t},
                               {'i': 0, 'j': 2, 'R': (0, -1), 't': t}])
    '''
    error_handling.positive_real(a, 'a')
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)},
                       {'tag': 'b', 'r0': (a, 0.)},
                       {'tag': 'c', 'r0': (0., a)}]
    return _lat(unit_cell, [(2*a, 0.), (0., 2*a)])
