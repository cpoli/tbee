"""
Density of states from a set of eigenenergies, real-space
(:class:`tbee.system.System`) or reciprocal-space
(:class:`tbee.kspace.KSpace`, sampled over a k-mesh).
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

import tbee.error_handling as error_handling


def density_of_states(
    energies: ArrayLike,
    e_grid: ArrayLike | None = None,
    broadening: float = 0.05,
    kernel: str = 'gaussian',
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    r'''
    Get the density of states, broadened by a Gaussian or Lorentzian
    kernel of width *broadening*:

    .. math::

        \rho(E) = \sum_n g(E-E_n)\, ,\quad
        g(x) = \frac{1}{\sqrt{2\pi}\sigma}e^{-x^2/2\sigma^2}\ \text{(gaussian)}
        \ \text{or}\
        g(x) = \frac{1}{\pi}\frac{\sigma}{x^2+\sigma^2}\ \text{(lorentzian)}

    Each level contributes a kernel of unit area, so
    :math:`\int\rho(E)dE` equals the number of levels in *energies*,
    for an *e_grid* wide enough to contain the tails.

    :param energies: Array of (real) eigenenergies. Any shape (e.g. the
        *en* attribute of **System**, or of **KSpace** after *get_bands*
        over a k-mesh -- flattened automatically).
    :param e_grid: Real ndarray. Default value None. Energies at which to
        evaluate the density of states. If None, a grid of 401 points
        spanning ``[min(energies)-3*broadening, max(energies)+3*broadening]``
        is used.
    :param broadening: Positive real number. Default value 0.05. Kernel width
        :math:`\sigma`.
    :param kernel: String. Default value 'gaussian'. 'gaussian' or 'lorentzian'.

    :returns:
        * **e_grid** -- Real ndarray. The energy grid used.
        * **dos** -- Real ndarray, same shape as *e_grid*. Density of states.
    '''
    error_handling.ndarray_empty(np.asarray(energies), 'energies')
    error_handling.positive_real(broadening, 'broadening')
    error_handling.dos_kernel(kernel)
    energies = np.asarray(energies).real.astype('f8').ravel()
    if e_grid is None:
        pad = 3 * broadening
        e_grid = np.linspace(energies.min() - pad, energies.max() + pad, 401)
    else:
        error_handling.ndarray_empty(np.asarray(e_grid), 'e_grid')
        e_grid = np.asarray(e_grid, dtype='f8')
    diff = e_grid[:, None] - energies[None, :]
    if kernel == 'gaussian':
        weight = np.exp(-diff**2 / (2*broadening**2)) / (broadening*np.sqrt(2*np.pi))
    else:
        weight = (broadening/np.pi) / (diff**2 + broadening**2)
    return e_grid, weight.sum(axis=1)
