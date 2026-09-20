r"""
A Topological Flat Band on the Kagome Lattice
===================================================

The kagome lattice's flat band (see :doc:`/api/gallery/flat_bands/plot_flat_bands`) isn't
isolated: at zero field it touches the middle dispersive band exactly
at :math:`\Gamma` (a symmetry-protected degeneracy, both at
:math:`E=-2t`). Making the nearest-neighbor hopping complex -- the same
phase :math:`t e^{i\phi}` on *every* bond, regardless of which of the
lattice's two triangle orientations it belongs to -- breaks the mirror
symmetries protecting that degeneracy (while leaving the threefold
rotation intact) and opens a genuine gap across the whole Brillouin
zone. The resulting lower band is no longer exactly flat, but it stays
topologically nontrivial: a Chern insulator with no net magnetic field,
built entirely from a real, physically motivated mechanism -- the
scalar spin chirality of a canted magnetic texture on the kagome
lattice acts, for the itinerant electrons, exactly like this complex
hopping (Ohgushi, Murakami, and Nagaosa, 2000).
"""
import numpy as np
import matplotlib.pyplot as plt

import tbee.lattices as lattices
from tbee.kspace import KSpace, reciprocal_vectors


lat = lattices.kagome()
t1 = 1.
b1, b2 = (np.array(v) for v in reciprocal_vectors(lat.prim_vec))
Gamma, K, M_pt = np.zeros(2), (b1 - b2) / 3, b1 / 2


def kagome_chiral(phi):
    '''Kagome lattice, nearest-neighbor hopping t1*exp(i*phi) on every bond.'''
    kag = KSpace(lat)
    t = t1 * np.exp(1j * phi)
    kag.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': t},
                            {'i': 0, 'j': 2, 'R': (0, 0), 't': t},
                            {'i': 0, 'j': 2, 'R': (0, -1), 't': t},
                            {'i': 1, 'j': 2, 'R': (0, 0), 't': t},
                            {'i': 1, 'j': 2, 'R': (1, -1), 't': t}])
    return kag


# %%
# phi=0: the flat band touches the middle band at Gamma
# ------------------------------------------------------------------

kag0 = kagome_chiral(phi=0.)
en_gamma_0 = np.sort(np.linalg.eigvalsh(kag0.get_ham(Gamma)))
print('phi=0, E(Gamma) = {} (bottom two exactly degenerate at -2t).'.format(np.round(en_gamma_0, 6)))
assert np.isclose(en_gamma_0[0], en_gamma_0[1], atol=1e-10)

# %%
# phi != 0: a real gap opens across the whole Brillouin zone
# ------------------------------------------------------------------

phi = 0.3
kag = kagome_chiral(phi)

nk = 60
ks = [(i / nk) * b1 + (j / nk) * b2 for i in range(nk) for j in range(nk)]
en_mesh = np.sort(np.array([np.linalg.eigvalsh(kag.get_ham(k)) for k in ks]), axis=1)
gap01 = en_mesh[:, 1].min() - en_mesh[:, 0].max()
gap12 = en_mesh[:, 2].min() - en_mesh[:, 1].max()
bandwidth0 = en_mesh[:, 0].max() - en_mesh[:, 0].min()
print('phi={}: gap below band 0 -> band 1 = {:.4f}, band 1 -> band 2 = {:.4f}, '
           'band 0 bandwidth = {:.4f}.'.format(phi, gap01, gap12, bandwidth0))
assert gap01 > 0.3
assert gap12 > 0.3

# %%
# Chern numbers: -1, 0, +1
# ------------------------------------------------------------------
# Exact flatness and a nonzero Chern number cannot coexist without
# further fine-tuning (this band's flatness ratio, gap/bandwidth, is
# only about 1 here -- the 2011 papers that proposed *nearly*-flat
# Chern bands as a lattice route to fractional Chern insulators pushed
# this ratio far higher via additional further-neighbor terms). What
# does survive intact is the topology.

chern = [kag.chern_number(bands=[n], nk=60) for n in range(3)]
print('Chern numbers (bottom, middle, top): {}'.format(np.round(chern, 4)))
assert np.isclose(chern[0], -1., atol=1e-2)
assert np.isclose(chern[1], 0., atol=1e-2)
assert np.isclose(chern[2], 1., atol=1e-2)
assert np.isclose(sum(chern), 0., atol=1e-2)
print('Bottom band is a C=-1 Chern insulator; the three bands sum to C=0. OK')

# %%
# Berry curvature of the isolated lower band
# ------------------------------------------------------------------

curv = kag.berry_curvature(bands=[0], nk=60)
fig, ax = plt.subplots()
im = ax.imshow(curv.T, origin='lower', extent=[0, 1, 0, 1], aspect='auto', cmap='RdBu')
ax.set_xlabel('$k_1$ (fractional)')
ax.set_ylabel('$k_2$ (fractional)')
ax.set_title('Berry curvature of the lower (Chern) band')
fig.colorbar(im, ax=ax)

# %%
# Band structure: flat-and-touching vs. gapped-and-topological
# --------------------------------------------------------------------

fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
for axis, model, title in [(ax1, kag0, r'$\phi=0$: flat band touches band 1'),
                                        (ax2, kag, r'$\phi={:.1f}$: gapped, $C=-1$'.format(phi))]:
    ks_dist, en = model.k_path([Gamma, K, M_pt, Gamma], nk=60)
    for n in range(3):
        axis.plot(ks_dist, en[:, n], 'b')
    axis.set_title(title)
    axis.set_xlabel('$k$')
ax1.set_ylabel('$E$')
fig2.tight_layout()
