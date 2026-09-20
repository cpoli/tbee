r"""
Hofstadter's Butterfly
===========================

Sweeping the flux per plaquette continuously (rather than fixing it, as
the Aharonov-Bohm ring example does) reveals a self-similar, fractal
spectrum: at rational flux p/q, each Landau-level-like band splits into q
sub-bands, and the whole pattern is periodic in the flux with period one
flux quantum. This script sweeps
:meth:`~tbee.system.System.set_magnetic_field` over a finite square-lattice
flake and plots the resulting "butterfly".
"""
import numpy as np
import matplotlib.pyplot as plt

from tbee.lattice import Lattice
from tbee.system import System


N1, N2 = 16, 16  # flake size


def flake(alpha):
    lat = Lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.), (0., 1.)])
    lat.get_lattice(n1=N1, n2=N2)
    sys = System(lat)
    sys.set_hopping([{'n': 1, 't': 1.}])
    if alpha != 0.:
        sys.set_magnetic_field(alpha=alpha)
    sys.get_ham()
    sys.get_eig()
    return sys


# %%
# Sanity checks that hold at every flux
# --------------------------------------------
# These hold regardless of the (fractal, hard to predict in detail)
# spectrum shape itself:
#
# 1. :math:`|E|\le 4|t|` always (Gershgorin: at most 4 neighbors of hopping
#    magnitude :math:`|t|=1`, and the Peierls substitution changes their
#    phase, never their magnitude).
# 2. The spectrum is symmetric under :math:`E \to -E` (the square lattice
#    is bipartite, and that chiral symmetry survives the Peierls
#    substitution since it only ever multiplies hoppings by a phase).
# 3. The spectrum is periodic in :math:`\alpha` with period 1 (one flux
#    quantum per plaquette is gauge-equivalent to no flux at all).

alphas = np.linspace(0., 1., 161)
spectra = [flake(alpha).en.real for alpha in alphas]

max_abs_e = max(np.max(np.abs(en)) for en in spectra)
print('max|E| over the whole sweep: {:.4f} (Gershgorin bound: 4.0000)'.format(max_abs_e))
assert max_abs_e <= 4. + 1e-8

sym_ok = all(np.allclose(np.sort(en), -np.sort(en)[::-1], atol=1e-6) for en in spectra)
print('Particle-hole symmetry E -> -E holds at every flux: {}'.format(sym_ok))
assert sym_ok

periodic_ok = np.allclose(np.sort(spectra[0]), np.sort(spectra[-1]), atol=1e-6)
print('Spectrum at alpha=0 matches alpha=1 (period-1 flux quantum): {}'.format(periodic_ok))
assert periodic_ok

# %%
# The butterfly itself
# ---------------------------

fig, ax = plt.subplots(figsize=(6, 7))
for alpha, en in zip(alphas, spectra):
    ax.plot(alpha*np.ones_like(en), en, '.', ms=0.5, color='b')
ax.set_xlabel(r'flux per plaquette $\alpha$ (in units of $\Phi_0$)')
ax.set_ylabel('$E$')
ax.set_title("Hofstadter's butterfly ({}x{} flake)".format(N1, N2))
