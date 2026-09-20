Tutorial
========

This walks through **tbee**'s pieces in the order you would normally use
them: define a lattice, build a real-space Tight-Binding model and
diagonalize it, build the same model's Bloch Hamiltonian in reciprocal
space, then layer on more physics (disorder, strain, a magnetic field,
spin-orbit coupling, topology, edge states). Runnable versions of most of
this live in the ``examples/`` directory, referenced throughout.

A Tight-Binding Hamiltonian describes electrons hopping between localized
orbitals:

.. math::

    H = \sum_i \epsilon_i c_i^\dagger c_i + \sum_{i \neq j} t_{ij}\, c_i^\dagger c_j

where :math:`\epsilon_i` are onsite energies and :math:`t_{ij}` are
hopping amplitudes between orbitals *i* and *j*. Everything in **tbee**
is about building this matrix (or, in reciprocal space, its Bloch analog
:math:`H(\mathbf{k})`) and then doing something with it: diagonalizing
it, plotting it, evolving a wavepacket under it.


Building a lattice
-------------------

:class:`tbee.lattice.Lattice` defines the geometry: a unit cell (a list of
orbitals, each with a sublattice *tag* and a position) and one or two
primitive vectors::

    from tbee.lattice import Lattice

    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}]
    prim_vec = [(1., 0.), (0., 1.)]
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=6, n2=6)   # 6x6 unit cells -> 36 sites

``lat.coor`` now holds every site's position and tag. Sculpt the shape
with methods like ``remove_sites``, ``ellipse_in``/``ellipse_out``,
``boundary_line``, or by adding/subtracting two lattices (``lat1 + lat2``).
:mod:`tbee.lattices` has a handful of common lattices (chain, square,
triangular, honeycomb, kagome, Lieb) ready to use instead of writing out
*unit_cell*/*prim_vec* by hand.

See ``lat.plot()`` to look at what you built, and
``examples/examples_lattice.ipynb`` for many more shapes.


Real space: build and solve a finite system
---------------------------------------------

:class:`tbee.system.System` takes a *Lattice* and builds the Hamiltonian
from a list of hoppings::

    from tbee.system import System

    sys = System(lat)
    sys.set_onsite({'a': 0.})
    sys.set_hopping([{'n': 1, 't': 1.}])   # nearest-neighbor hopping t=1
    sys.get_ham()
    sys.get_eig()
    print(sys.en)   # eigenenergies, sorted

Hoppings are addressed by neighbor order (``'n': 1`` is nearest-neighbor,
``'n': 2`` next-nearest, ...; ``sys.print_distances()`` lists them), and
optionally further filtered by bond angle or by sublattice-pair tag --
see :meth:`tbee.system.System.set_hopping`'s docstring for the full
mini-language. For anything that doesn't fit that pattern, ``i,j``-indexed
hoppings work directly on the actual real-space site indices printed by
``lat.plot(plt_index=True)``. Hamiltonians can be Hermitian or, if you fill
in both the *i<j* and *i>j* parts independently, non-Hermitian.

``tbee.plot.Plot`` plots the lattice, the spectrum, sublattice
polarization, eigenstate intensities, and the density of states
(:meth:`~tbee.plot.Plot.dos`)::

    from tbee.plot import Plot
    p = Plot(sys)
    p.spectrum()
    p.dos(broadening=0.1)

.. minigallery:: ../../examples/tight_binding/plot_visualizing_a_model.py


Reciprocal space: Bloch Hamiltonians and band structures
------------------------------------------------------------

A finite flake is one way to look at a lattice; the other is to keep it
infinite and periodic, and work with its Bloch Hamiltonian
:math:`H(\mathbf{k})`. :class:`tbee.kspace.KSpace` builds it from a small
set of *intra-unit-cell* hoppings, each tagged by which neighboring cell
(:math:`\mathbf{R} = n_1\mathbf{a}_1+n_2\mathbf{a}_2`) it connects to::

    from tbee.kspace import KSpace, reciprocal_vectors
    import tbee.lattices as lattices
    import numpy as np

    lat = lattices.honeycomb()
    gra = KSpace(lat)
    gra.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}])

Only one representative of each bond is needed -- the reverse bond
(:math:`j\to i`, :math:`\mathbf{R}\to-\mathbf{R}`) is added automatically
as its Hermitian conjugate. Diagonalize at a single k-point with
``get_ham``/``LA.eigh``, over an explicit set of k-points with
``get_bands``, or along a path through high-symmetry points with
``k_path`` (:func:`~tbee.kspace.reciprocal_vectors` gives the reciprocal
lattice vectors to build that path from)::

    b1, b2 = (np.array(v) for v in reciprocal_vectors(lat.prim_vec))
    Gamma, K, M = np.zeros(2), (b1 - b2)/3, b1/2
    gra.k_path([Gamma, K, M, Gamma], nk=60)
    fig = gra.plot_bands(node_labels=[r'$\Gamma$', 'K', 'M', r'$\Gamma$'])

``gra.mesh_bands(nk)`` samples a uniform mesh over the whole Brillouin
zone instead (used internally by ``plot_dos``, and for the topology tools
below). See ``examples/tight_binding/plot_graphene_bands.py``.

.. minigallery:: ../../examples/tight_binding/plot_graphene_bands.py


Disorder, strain, and a magnetic field
-----------------------------------------

On a *System*, ``set_onsite_dis``/``set_hopping_dis`` add uniform random
disorder; ``GrapheneSystem.set_hop_linear_strain`` (in :mod:`tbee.graphene`)
applies triaxial strain. An orbital magnetic field is added via the
Peierls substitution -- each hopping amplitude gets multiplied by a phase
equal to the line integral of the vector potential along the bond::

    sys.set_hopping_manual(hop_dict)
    sys.set_magnetic_field(alpha=0.01)   # alpha = B/Phi_0, flux quanta per unit area
    sys.get_ham()

``set_peierls_phase`` takes an arbitrary phase function for a
non-uniform field or a different gauge. See
``examples/magnetic_field/plot_magnetic_field.py`` for an Aharonov-Bohm ring, where the
spectrum is checked to be exactly periodic in the enclosed flux.

.. minigallery:: ../../examples/magnetic_field/plot_magnetic_field.py


Spin and spin-orbit coupling
--------------------------------

``KSpace(lat, spin=True)`` gives every site a spin-1/2 degree of freedom
(orbitals double: ``2*i, 2*i+1`` are the up/down components of site *i*).
``set_onsite``/``set_hopping`` then also accept 2x2 matrices -- built from
:data:`tbee.kspace.PAULI`'s Pauli matrices -- for spin-dependent terms::

    from tbee.kspace import PAULI

    kmele = KSpace(lat, spin=True)
    kmele.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                              {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                              {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}])
    for R in [(1, 0), (0, 1), (1, -1)]:
        kmele.set_hopping([{'i': 0, 'j': 0, 'R': R, 't': 1j*lam*PAULI['z']}])   # intrinsic SOC
        kmele.set_hopping([{'i': 1, 'j': 1, 'R': R, 't': -1j*lam*PAULI['z']}])

A plain number is still accepted for a spin-independent hopping (expanded
to ``t*PAULI['0']``), and onsite values accept a pair ``(E_up, E_down)``
for a Zeeman-like splitting.


Topology: Berry curvature and Chern numbers
------------------------------------------------

For a periodic (*KSpace*) model, :meth:`~tbee.kspace.KSpace.chern_number`
gives the Chern number of a group of bands (an isolated band, or several
occupied bands below a gap), computed by the gauge-invariant
Fukui-Hatsugai-Suzuki lattice method -- an integer, for a group of bands
with a gap above and below it everywhere in the Brillouin zone::

    chern = gra.chern_number(bands=[0], nk=40)    # ~0 for plain graphene (gapless!)

:meth:`~tbee.kspace.KSpace.berry_curvature` gives the underlying
:math:`\mathbf{k}`-resolved curvature (``chern_number`` is just its sum,
divided by :math:`2\pi`). See ``examples/topology/plot_haldane_topology.py`` for the
Haldane model -- the original Chern insulator -- reproducing its
topological phase transition (Chern number 1 in one phase, 0 in the
other) as a sublattice mass term is tuned.

.. minigallery:: ../../examples/topology/plot_haldane_topology.py


Edge states: cutting a ribbon
----------------------------------

Bulk topology has a physical consequence at an edge. :func:`tbee.kspace.ribbon`
cuts a ribbon -- periodic along one primitive vector, finite (open
boundary, a chosen number of unit cells) along the other -- out of any
periodic model, keeping its hoppings intact::

    from tbee.kspace import ribbon

    list_hop = [{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                      {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                      {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}]
    rib = ribbon(lat, list_hop, width=30, direction=1)
    rib.k_path([(-np.pi,), (np.pi,)], nk=200)
    fig = rib.plot_bands()

For a zigzag graphene ribbon this reproduces the famous zero-energy edge
flat band; cutting a ribbon from the Kane-Mele model above instead shows
helical edge states crossing the bulk gap, protected by Kramers' theorem.
See ``examples/topology/plot_edge_states.py``, where both are checked numerically (the
second also against a large real-space flake, as an independent
cross-check).

.. minigallery:: ../../examples/topology/plot_edge_states.py


Where next
--------------

* :doc:`tbee` -- the full API reference, generated from the source
  docstrings.
* ``examples/`` in the repository -- runnable scripts and notebooks for
  everything above, plus more (kagome, Lieb, propagation, ...).
* ``tests/`` -- every claim above (flat bands, Chern numbers, Kramers
  degeneracy, flux periodicity, ...) is checked against an analytic or
  independently-computed result somewhere in the test suite; reading them
  alongside the code they test is a good way to see the underlying physics
  made concrete.
