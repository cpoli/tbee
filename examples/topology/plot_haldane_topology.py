r"""
The Haldane Model: Berry Curvature and a Topological Phase Transition
===========================================================================

The Haldane model is honeycomb graphene with a complex next-nearest-neighbor
hopping :math:`t_2 e^{i\phi}` (breaking time-reversal symmetry, e.g. via a
staggered flux pattern with zero net flux) and a staggered sublattice
onsite energy :math:`\pm M` (breaking inversion symmetry). It is the first
model shown to realize a Chern insulator -- a gapped phase with quantized
Hall conductance and no net magnetic field.

The topological/trivial phase boundary sits at
:math:`|M| = \sqrt{3}\,t_2\,|\sin\phi|` (for this particular choice of
which 3 next-nearest-neighbor vectors carry the phase :math:`+\phi` vs.
:math:`-\phi` -- the prefactor is convention-dependent and was pinned down
numerically below, rather than assumed): the lower band's Chern number,
computed via :meth:`~tbee.kspace.KSpace.chern_number`, is
:math:`\pm 1` for :math:`|M|` below that, 0 above it.
"""
import numpy as np
import matplotlib.pyplot as plt

from tbee.lattice import Lattice
from tbee.kspace import KSpace, reciprocal_vectors


DX, DY = 0.5 * 3 ** 0.5, 0.5
unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
prim_vec = [(2*DX, 0.), (DX, 1.5)]
t1, t2, phi = 1., 0.2, np.pi / 2
M_c = np.sqrt(3) * t2 * abs(np.sin(phi))  # critical mass (see module docstring)


def haldane(M):
    '''Build the Haldane model with staggered onsite energy +-M.'''
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    hal = KSpace(lat)
    hal.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': t1}])
    # next-nearest-neighbor hopping: same 3 lattice vectors for both
    # sublattices, but with opposite chirality (t2*exp(+-i*phi)) -- this
    # circulating "staggered flux" is what breaks time-reversal symmetry
    # without any net magnetic field through the unit cell.
    for R in [(1, 0), (0, 1), (1, -1)]:
        hal.set_hopping([{'i': 0, 'j': 0, 'R': R, 't': t2*np.exp(1j*phi)}])
        hal.set_hopping([{'i': 1, 'j': 1, 'R': R, 't': t2*np.exp(-1j*phi)}])
    hal.set_onsite({'a': M, 'b': -M})
    return hal


# %%
# Chern number across the topological phase transition
# ------------------------------------------------------------

masses = np.linspace(0., 2*M_c, 21)
chern = [haldane(M).chern_number(bands=[0], nk=40) for M in masses]

print('Critical mass M_c = sqrt(3)*t2*sin(phi) = {:.4f}'.format(M_c))
print('Chern number at M=0            (topological): {:.4f}'.format(chern[0]))
print('Chern number at M=2*M_c        (trivial):      {:.4f}'.format(chern[-1]))
assert np.isclose(chern[0], 1., atol=1e-2)
assert np.isclose(chern[-1], 0., atol=1e-2)
print('Phase transition reproduced: C = 1 (topological) -> C = 0 (trivial). OK')

fig, ax = plt.subplots()
ax.plot(masses/M_c, chern, 'o-b')
ax.axvline(1., color='k', ls='--', lw=1)
ax.set_xlabel('$M/M_c$')
ax.set_ylabel('Chern number (lower band)')
ax.set_title('Haldane model: topological phase transition')

# %%
# Berry curvature in the topological phase
# ------------------------------------------------
# :meth:`~tbee.kspace.KSpace.berry_curvature` concentrates near the Dirac
# points (where the gap is smallest), with total flux :math:`2\pi`.

hal_topological = haldane(M=0.)
curv = hal_topological.berry_curvature(bands=[0], nk=60)
fig2, ax2 = plt.subplots()
im = ax2.imshow(curv.T, origin='lower', extent=[0, 1, 0, 1], aspect='auto', cmap='RdBu')
ax2.set_xlabel('$k_1$ (fractional)')
ax2.set_ylabel('$k_2$ (fractional)')
ax2.set_title('Berry curvature of the lower band')
fig2.colorbar(im, ax=ax2)

# %%
# Band structure in the topological phase
# ------------------------------------------------

b1, b2 = (np.array(v) for v in reciprocal_vectors(prim_vec))
Gamma, K, M_pt = np.zeros(2), (b1 - b2) / 3, b1 / 2
hal_topological.k_path([Gamma, K, M_pt, Gamma], nk=60)
fig3 = hal_topological.plot_bands(node_labels=[r'$\Gamma$', 'K', 'M', r'$\Gamma$'])
