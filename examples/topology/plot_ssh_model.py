r"""
The SSH Model: Bulk-Boundary Correspondence
=================================================

A 1D dimerized chain -- alternating intracell hopping v and intercell
hopping w -- is the earliest and simplest example of a topological
insulator: for w > v the chain is topological and an *open* finite chain
hosts a pair of exponentially localized, near-zero-energy edge modes, one
per end. For v > w it is trivial and no such modes exist. The bulk gap
(:math:`2|v-w|`) closes exactly at the transition, v = w.
"""
import numpy as np
import matplotlib.pyplot as plt

from tbee.lattice import Lattice
from tbee.system import System
from tbee.kspace import KSpace


unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (0.5, 0.)}]
prim_vec = [(1., 0.)]


# %%
# Bulk (k-space): the gap closes at v = w
# --------------------------------------------
# The two bands are :math:`E(k)=\pm|v+we^{-ik}|`, so the gap is
# :math:`2\min_k|v+we^{-ik}| = 2|v-w|` (minimized at :math:`k=\pi`),
# closing exactly at the transition v = w.

def ssh_bulk(v, w):
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    ssh = KSpace(lat)
    ssh.set_hopping([{'i': 0, 'j': 1, 'R': (0,), 't': v},
                            {'i': 0, 'j': 1, 'R': (-1,), 't': w}])
    return ssh


for v, w in [(1., 0.5), (0.5, 1.), (1., 1.)]:
    en = ssh_bulk(v, w).mesh_bands(nk=200)
    gap = en[:, 1].min() - en[:, 0].max()
    assert np.isclose(gap, 2*abs(v - w), atol=1e-6)
print('Bulk gap = 2|v - w| confirmed for several (v, w).')

ssh = ssh_bulk(v=0.6, w=1.)
ssh.k_path([(-np.pi,), (np.pi,)], nk=200)
fig = ssh.plot_bands()

# %%
# Open finite chain: edge states in the topological phase only
# --------------------------------------------------------------------
# Edge states appear only for w > v, decaying exponentially into the bulk
# from each end.

def ssh_open_chain(v, w, n_cells):
    lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    lat.get_lattice(n1=n_cells)
    sys = System(lat)
    hop = {}
    for n in range(n_cells):
        hop[(2*n, 2*n + 1)] = v
        if n < n_cells - 1:
            hop[(2*n + 1, 2*n + 2)] = w
    sys.set_hopping_manual(hop)
    sys.get_ham()
    return sys


n_cells = 30
sys_trivial = ssh_open_chain(v=1., w=0.5, n_cells=n_cells)
sys_trivial.get_eig()
n_zero_trivial = np.sum(np.abs(sys_trivial.en.real) < 1e-6)
print('Trivial phase (v>w): {} exactly-zero-energy states (expect 0).'.format(n_zero_trivial))
assert n_zero_trivial == 0

sys_topological = ssh_open_chain(v=0.5, w=1., n_cells=n_cells)
sys_topological.get_eig(eigenvec=True)
n_zero_topological = np.sum(np.abs(sys_topological.en.real) < 1e-6)
print('Topological phase (w>v): {} exactly-zero-energy states (expect 2).'.format(n_zero_topological))
assert n_zero_topological == 2

# The two near-zero modes must be exponentially localized at the two ends.
idx = np.argsort(np.abs(sys_topological.en.real))[:2]
intensity = np.abs(sys_topological.rn[:, idx[0]]) ** 2
left_weight = intensity[:6].sum()
bulk_weight = intensity[2*n_cells//2 - 3: 2*n_cells//2 + 3].sum()
print('Edge-mode weight on the first 6 sites: {:.4f}; on 6 sites at the '
           'chain center: {:.6f}.'.format(left_weight, bulk_weight))
assert left_weight > 100 * bulk_weight

fig2, axes = plt.subplots(1, 2, figsize=(10, 4))
for ax, sys_, title in zip(axes, [sys_trivial, sys_topological],
                                       ['Trivial (v > w)', 'Topological (w > v)']):
    ax.plot(np.sort(sys_.en.real), 'o', ms=3)
    ax.axhline(0., color='k', lw=0.5)
    ax.set_title(title)
    ax.set_xlabel('state index')
axes[0].set_ylabel('$E$')
fig2.suptitle('SSH open chain ({} unit cells): spectrum'.format(n_cells))
