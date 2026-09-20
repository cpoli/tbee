r"""
The Thouless Quantum Pump (Rice-Mele Model)
=================================================

Take the SSH model (:doc:`plot_ssh_model`) and cycle its two parameters --
the dimerization v-w and a staggered onsite potential -- slowly around a
closed loop enclosing the SSH gap-closing point. Thouless showed that the
charge transported around one full cycle is exactly quantized, equal to
the Chern number of the :math:`(k,\phi)` torus swept out by the cycle -- a
1D adiabatic pump is, in this sense, a slice through a 2D Chern insulator,
with the pump parameter :math:`\phi` playing the role of a second crystal
momentum.

.. math::

    H(k, \phi) = \begin{pmatrix} \Delta\sin\phi & v(\phi)+w(\phi)e^{-ik}
    \\ \text{h.c.} & -\Delta\sin\phi \end{pmatrix}, \qquad
    v(\phi) = v_0 + \delta\cos\phi,\ \ w(\phi) = v_0 - \delta\cos\phi

Because :math:`v(\phi)` and :math:`w(\phi)` only ever appear as
:math:`\cos\phi = (e^{i\phi}+e^{-i\phi})/2`, this is a completely ordinary
2D Bloch Hamiltonian, built with :class:`tbee.kspace.KSpace` exactly like
any other -- no new machinery needed, just one extra "synthetic"
reciprocal dimension.
"""
import numpy as np
import matplotlib.pyplot as plt

from tbee.lattice import Lattice
from tbee.kspace import KSpace


unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (0.5, 0.)}]
prim_vec = [(1., 0.), (0., 1.)]  # 2nd direction: the synthetic pump parameter


def rice_mele(v0, delta, Delta, onsite_offset=0.):
    '''
    Rice-Mele model as a 2D Bloch Hamiltonian H(k, phi): a family of SSH
    chains (dimerization v(phi), w(phi)) with a staggered onsite potential
    Delta*sin(phi) + onsite_offset, phi playing the role of a second,
    synthetic crystal momentum.
    '''
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    rm = KSpace(lat)
    rm.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': v0},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': v0},
                            {'i': 0, 'j': 1, 'R': (0, 1), 't': delta/2},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': delta/2},
                            {'i': 0, 'j': 1, 'R': (-1, 1), 't': -delta/2},
                            {'i': 0, 'j': 1, 'R': (-1, -1), 't': -delta/2},
                            {'i': 0, 'j': 0, 'R': (0, 1), 't': Delta/(2j)},
                            {'i': 1, 'j': 1, 'R': (0, 1), 't': -Delta/(2j)}])
    if onsite_offset:
        rm.set_onsite({'a': onsite_offset, 'b': -onsite_offset})
    return rm


v0, delta, Delta = 1., 0.5, 0.6


# %%
# The pump is quantized only when the loop encircles the gap-closing point
# -----------------------------------------------------------------------------
# Chern number :math:`\pm 1` when the (v-w, staggering) loop encircles the
# SSH gap-closing point at the origin; offsetting the loop so it misses
# the origin entirely gives a trivial, unquantized (here, exactly zero)
# pump instead.

rm = rice_mele(v0, delta, Delta)
chern = rm.chern_number(bands=[0], nk=50)
print('Loop encircling the origin: Chern number = {:.4f} (expect exactly +-1).'.format(chern))
assert np.isclose(abs(chern), 1., atol=1e-3)

rm_offset = rice_mele(v0, delta, Delta, onsite_offset=2*Delta)
chern_offset = rm_offset.chern_number(bands=[0], nk=50)
print('Loop missing the origin (large constant offset): Chern number = {:.4f} '
           '(expect exactly 0).'.format(chern_offset))
assert np.isclose(chern_offset, 0., atol=1e-6)

# %%
# The physical signature: polarization winds by exactly one lattice vector
# ------------------------------------------------------------------------------
# The instantaneous electric polarization (the lower band's Zak/Berry
# phase at fixed phi, in units of the lattice constant) winds by exactly
# one full lattice vector over one pump cycle -- "one electron pumped
# through the bulk per cycle".

def polarization(rm, phi, nk=300):
    ks = np.linspace(0., 2*np.pi, nk, endpoint=False)
    us = [np.linalg.eigh(rm.get_ham((k, phi)))[1][:, 0] for k in ks]
    prod = np.prod([np.vdot(us[n], us[(n + 1) % nk]) for n in range(nk)])
    return -np.angle(prod) / (2*np.pi)


phis = np.linspace(0., 2*np.pi, 41)
pols = np.array([polarization(rm, phi) for phi in phis])
pols_unwrapped = np.unwrap(pols * 2*np.pi) / (2*np.pi)
winding = pols_unwrapped[-1] - pols_unwrapped[0]
print('Polarization winding over one pump cycle: {:.4f} (expect exactly +-1).'.format(winding))
assert np.isclose(abs(winding), 1., atol=1e-2)

fig, ax = plt.subplots()
ax.plot(phis / (2*np.pi), pols_unwrapped - pols_unwrapped[0], 'o-b')
ax.set_xlabel(r'pump cycle $\phi / 2\pi$')
ax.set_ylabel('polarization (lattice constants)')
ax.set_title('Thouless pump: polarization vs. pump parameter')
