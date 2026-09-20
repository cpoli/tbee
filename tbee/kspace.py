from __future__ import annotations

from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import scipy.linalg as LA
import tbee.error_handling as error_handling
import tbee.dos as dos
from tbee.lattice import Lattice


PI = np.pi

#: Pauli matrices (plus the identity, key ``'0'``), for building spinful
#: hoppings/onsite terms (spin-orbit coupling, Zeeman splitting, ...) when
#: :class:`KSpace` is constructed with ``spin=True``.
PAULI = {
    '0': np.eye(2, dtype='c16'),
    'x': np.array([[0., 1.], [1., 0.]], dtype='c16'),
    'y': np.array([[0., -1j], [1j, 0.]], dtype='c16'),
    'z': np.array([[1., 0.], [0., -1.]], dtype='c16'),
}


#################################
# CLASS KSPACE
#################################


def reciprocal_vectors(prim_vec: list[tuple[float, float]]) -> list[tuple[float, float]]:
    r'''
    Get the reciprocal lattice vectors :math:`\mathbf{b}_i` such that
    :math:`\mathbf{a}_i\cdot\mathbf{b}_j = 2\pi\delta_{ij}`.

    :param prim_vec: List of one/two tuples. Primitive vectors (see class **lattice**).

    :returns:
        * **rec_vec** -- List of one/two tuples. Reciprocal vectors.
    '''
    if len(prim_vec) == 1:
        ax, ay = prim_vec[0]
        norm2 = ax ** 2 + ay ** 2
        return [(2*PI*ax/norm2, 2*PI*ay/norm2)]
    (a1x, a1y), (a2x, a2y) = prim_vec
    area = a1x * a2y - a1y * a2x
    b1 = (2*PI*a2y/area, -2*PI*a2x/area)
    b2 = (-2*PI*a1y/area, 2*PI*a1x/area)
    return [b1, b2]


class KSpace():
    r'''
    Build and solve the Tight-Binding Bloch Hamiltonian :math:`H(\mathbf{k})`
    of a periodic lattice defined by the class **lattice**.

    Hoppings are defined between orbitals of the unit cell, separated by a
    lattice vector :math:`\mathbf{R} = n_1\mathbf{a}_1+n_2\mathbf{a}_2`:

    .. math::

        H_{ij}(\mathbf{k}) = \sum_{\mathbf{R}} t_{ij}(\mathbf{R})\,
        e^{i\mathbf{k}\cdot\mathbf{R}}

    :param lat: **lattice** class instance. Only *unit_cell* and *prim_vec*
        are used (the instance need not call *get_lattice*).
    :param spin: Boolean. Default value False. If True, every site of
        *unit_cell* carries a spin-1/2 degree of freedom (*norb* doubles to
        ``2*len(unit_cell)``, ordered site-major: orbitals ``2*i, 2*i+1``
        are the up/down components of site *i*). *set_onsite* and
        *set_hopping* then accept 2x2 (spin) matrices in addition to plain
        numbers, to build spin-orbit coupling or Zeeman terms -- see
        :data:`PAULI` for ready-made Pauli matrices.

    Example usage::

        # graphene, nearest-neighbor hopping t
        DX, DY = 0.5 * 3 ** 0.5, 0.5
        unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
        prim_vec = [(2*DX, 0.), (DX, 1.5)]
        lat = Lattice(unit_cell=unit_cell, prim_vec=prim_vec)
        gra = KSpace(lat)
        gra.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                                {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                                {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}])
    '''

    def __init__(self, lat: Lattice, spin: bool = False) -> None:
        error_handling.lat(lat)
        error_handling.boolean(spin, 'spin')
        self.lat = lat
        self.dim = len(lat.prim_vec)
        self.spin = spin
        self.n_sites = len(lat.unit_cell)
        self.norb = 2*self.n_sites if spin else self.n_sites
        self.tags = np.array([dic['tag'] for dic in lat.unit_cell])
        self.onsite = np.zeros(self.norb, 'c16')
        self._hop = []  # list of (i, j, R_cartesian (np.ndarray), t)
        self.rec_vec = reciprocal_vectors(lat.prim_vec)
        self.ks = np.array([])  # k-points of the last band-structure calculation
        self.ks_dist = np.array([])  # cumulative distance along the k-path
        self.nodes = np.array([])  # positions, along ks_dist, of the k-path nodes
        self.en = np.array([])  # bands, shape (len(ks), norb)

    def set_onsite(self, dict_onsite: dict[str, complex | Sequence[complex]]) -> None:
        '''
        Set the onsite energies, by sublattice tag.

        :param dict_onsite: Dictionary. key: tag, val: onsite energy
            (a plain number), or, if ``spin=True``, either a plain number
            (applied equally to both spins) or a pair ``(E_up, E_down)`` of
            numbers (a spin splitting, e.g. a Zeeman term along z).

        Example usage::

            kag.set_onsite({'a': 1., 'b': -1.})
            # spinful: same onsite energy for both spins on 'a', a Zeeman
            # splitting on 'b':
            kag_spin.set_onsite({'a': 1., 'b': (1., -1.)})
        '''
        error_handling.set_onsite_kspace(dict_onsite, self.lat.tags, self.spin)
        for tag, val in dict_onsite.items():
            sites = np.where(self.tags == tag)[0]
            if self.spin:
                e_up, e_down = (val, val) if isinstance(val, (int, float, complex)) else val
                self.onsite[2*sites] = e_up
                self.onsite[2*sites + 1] = e_down
            else:
                self.onsite[sites] = val

    def set_hopping(self, list_hop: list[dict]) -> None:
        r'''
        Set the hoppings between orbitals of the unit cell.

        Only one representative of each hopping needs to be given: its
        Hermitian conjugate (:math:`j\to i`, :math:`\mathbf{R}\to-\mathbf{R}`)
        is added automatically.

        :param list_hop: List of dictionaries with keys ('i', 'j', 'R', 't'):

            * 'i', 'j': Positive integers. Site indices within the unit cell
              (following the order of *unit_cell*).
            * 'R': Tuple of one/two integers :math:`(n_1, n_2)`. Lattice vector
              :math:`\mathbf{R}=n_1\mathbf{a}_1+n_2\mathbf{a}_2` separating the
              two sites.
            * 't': Complex number, or, if ``spin=True``, either a complex
              number (spin-independent hopping) or a 2x2 complex matrix (a
              general, possibly spin-mixing, hopping -- e.g. built from
              :data:`PAULI` for Rashba or intrinsic spin-orbit coupling).

        Example usage::

            # 1D chain, nearest-neighbor hopping t between the only orbital
            # and its right neighbor:
            chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
            # spinful: spin-independent hopping t, plus a Rashba-like
            # spin-flip term of strength alpha:
            chain_spin.set_hopping([{'i': 0, 'j': 0, 'R': (1,),
                                                    't': t*PAULI['0'] + 1j*alpha*PAULI['y']}])
        '''
        error_handling.set_hopping_kspace(list_hop, self.n_sites, self.dim, self.spin)
        for dic in list_hop:
            R_cart = np.zeros(2)
            for n, a in zip(dic['R'], self.lat.prim_vec):
                R_cart += n * np.array(a)
            i, j, t = dic['i'], dic['j'], dic['t']
            if self.spin:
                block = t*PAULI['0'] if isinstance(t, (int, float, complex)) else np.asarray(t, 'c16')
                for a in range(2):
                    for b in range(2):
                        self._hop.append((2*i+a, 2*j+b, R_cart, block[a, b]))
                        self._hop.append((2*j+b, 2*i+a, -R_cart, np.conj(block[a, b])))
            else:
                self._hop.append((i, j, R_cart, t))
                if not (i == j and not np.any(R_cart)):
                    self._hop.append((j, i, -R_cart, np.conj(t)))

    def clear_hopping(self) -> None:
        '''
        Clear the hoppings set by *set_hopping*.
        '''
        self._hop = []

    def get_ham(self, k: ArrayLike) -> NDArray[np.complex128]:
        r'''
        Get the dense Bloch Hamiltonian :math:`H(\mathbf{k})`.

        :param k: Tuple/list/ndarray of one/two real numbers. :math:`\mathbf{k}` point,
            in the same Cartesian frame as *prim_vec*.

        :returns:
            * **ham** -- Complex ndarray, shape (norb, norb).
        '''
        error_handling.k_vector(k, 'k', self.dim)
        k_cart = np.zeros(2)
        k_cart[:self.dim] = k
        ham = np.diag(self.onsite).astype('c16')
        for i, j, R_cart, t in self._hop:
            ham[i, j] += t * np.exp(1j * np.dot(k_cart, R_cart))
        return ham

    def get_bands(
        self, ks: ArrayLike, eigenvec: bool = False,
    ) -> NDArray[np.float64] | tuple[NDArray[np.float64], NDArray[np.complex128]]:
        r'''
        Diagonalize :math:`H(\mathbf{k})` over a set of k-points.

        :param ks: ndarray, shape (nk, dim). k-points.
        :param eigenvec: Boolean. Default value False. If True, also return
            the eigenvectors.

        :returns:
            * **en** -- Real ndarray, shape (nk, norb). Band energies, sorted ascending.
            * **vn** -- Complex ndarray, shape (nk, norb, norb), only if *eigenvec* is True.
              vn[k, :, n] is the nth eigenvector at ks[k].
        '''
        ks = np.atleast_2d(np.asarray(ks, dtype='f8'))
        self.ks = ks
        self.en = np.zeros((len(ks), self.norb))
        if eigenvec:
            vn = np.zeros((len(ks), self.norb, self.norb), 'c16')
        for i, k in enumerate(ks):
            ham = self.get_ham(k)
            if eigenvec:
                en, v = LA.eigh(ham)
                vn[i] = v
            else:
                en = LA.eigvalsh(ham)
            self.en[i] = en
        if eigenvec:
            return self.en, vn
        return self.en

    def k_path(
        self, points: list[ArrayLike], nk: int,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        r'''
        Build a k-path through a list of high-symmetry points, and get the
        associated bands.

        :param points: List of at least two k-points (each a tuple/list of
            one/two real numbers).
        :param nk: Positive integer. Number of k-points per path segment.

        :returns:
            * **ks_dist** -- Real ndarray. Cumulative distance along the path,
              to be used as the x-axis of a band-structure plot.
            * **en** -- Real ndarray, shape (len(ks_dist), norb). Band energies.
        '''
        error_handling.k_path_points(points, self.dim)
        error_handling.positive_int(nk, 'nk')
        points = np.atleast_2d(np.asarray(points, dtype='f8'))
        segments = [np.linspace(points[i], points[i+1], nk, endpoint=False)
                          for i in range(len(points) - 1)]
        ks = np.concatenate(segments + [points[-1:]])
        steps = np.linalg.norm(np.diff(ks, axis=0), axis=1)
        self.ks_dist = np.concatenate([[0.], np.cumsum(steps)])
        self.nodes = self.ks_dist[::nk][:len(points)-1].tolist() + [self.ks_dist[-1]]
        en = self.get_bands(ks)
        return self.ks_dist, en

    def mesh_grid(
        self, nk: int | tuple[int, int],
    ) -> tuple[list[NDArray[np.float64]], NDArray[np.float64]]:
        '''
        Private method. Build a uniform grid of fractional coordinates
        spanning the Brillouin zone (each in [0, 1)), and the corresponding
        Cartesian k-points.

        :param nk: Positive integer, or tuple of *dim* positive integers.
            Number of k-points along each reciprocal lattice vector.

        :returns:
            * **fracs** -- List of *dim* real ndarrays, shape (nk1, nk2) each
              (or (nk1,) in 1D): fractional coordinates of the grid.
            * **ks** -- Real ndarray, shape (nk1*nk2, dim) (or (nk1, dim) in 1D).
        '''
        error_handling.nk(nk, self.dim)
        if isinstance(nk, int):
            nk = (nk,) * self.dim
        rec_vec = [np.array(b) for b in self.rec_vec]
        if self.dim == 1:
            f1 = np.arange(nk[0]) / nk[0]
            ks = f1[:, None] * rec_vec[0][None, :self.dim]
            return [f1], ks
        f1, f2 = np.meshgrid(np.arange(nk[0])/nk[0], np.arange(nk[1])/nk[1], indexing='ij')
        ks = (f1.ravel()[:, None] * rec_vec[0][None, :]
                    + f2.ravel()[:, None] * rec_vec[1][None, :])
        return [f1, f2], ks

    def mesh_bands(self, nk: int | tuple[int, int]) -> NDArray[np.float64]:
        '''
        Diagonalize :math:`H(\\mathbf{k})` over a uniform mesh spanning the
        Brillouin zone.

        :param nk: Positive integer, or tuple of *dim* positive integers.
            Number of k-points along each reciprocal lattice vector.

        :returns:
            * **en** -- Real ndarray, shape (nk1*nk2, norb). Band energies
              over the mesh (flattened).
        '''
        _, ks = self.mesh_grid(nk)
        return self.get_bands(ks)

    def berry_curvature(
        self, bands: int | list[int], nk: int | tuple[int, int] = 30,
    ) -> NDArray[np.float64]:
        r'''
        Get the Berry curvature of a group of bands over a uniform
        Brillouin-zone mesh, using the gauge-invariant lattice method of
        Fukui, Hatsugai and Suzuki (J. Phys. Soc. Jpn. 74, 1674 (2005)):
        the flux through each mesh plaquette is minus the phase of the
        product of the (Slater-determinant) overlaps between the occupied
        subspaces at its four corners.

        :param bands: Positive integer, or list of positive integers. Band
            index, or indices of a group of bands (e.g. all occupied bands
            below a gap).
        :param nk: Positive integer, or tuple of 2 positive integers.
            Default value 30. Number of k-points along each reciprocal
            lattice vector.

        :returns:
            * **curv** -- Real ndarray, shape (nk1, nk2). Berry curvature
              (flux through each plaquette, in radians). Summing *curv* and
              dividing by :math:`2\pi` gives the Chern number, see
              *chern_number*.
        '''
        error_handling.dim_2(self.dim)
        if isinstance(bands, int):
            bands = [bands]
        error_handling.band_indices(bands, self.norb)
        if isinstance(nk, int):
            nk = (nk, nk)
        error_handling.nk(nk, 2)
        n1, n2 = nk
        _, ks = self.mesh_grid(nk)
        ks = ks.reshape(n1, n2, 2)
        v = np.zeros((n1, n2, self.norb, len(bands)), 'c16')
        for i1 in range(n1):
            for i2 in range(n2):
                _, vn = LA.eigh(self.get_ham(ks[i1, i2]))
                v[i1, i2] = vn[:, bands]
        curv = np.zeros((n1, n2))
        for i1 in range(n1):
            for i2 in range(n2):
                v1 = v[i1, i2]
                v2 = v[(i1+1) % n1, i2]
                v3 = v[(i1+1) % n1, (i2+1) % n2]
                v4 = v[i1, (i2+1) % n2]
                link = (np.linalg.det(v1.conj().T @ v2)
                             * np.linalg.det(v2.conj().T @ v3)
                             * np.linalg.det(v3.conj().T @ v4)
                             * np.linalg.det(v4.conj().T @ v1))
                curv[i1, i2] = -np.angle(link)
        return curv

    def chern_number(self, bands: int | list[int], nk: int | tuple[int, int] = 30) -> float:
        r'''
        Get the Chern number of a group of bands:

        .. math::

            C = \frac{1}{2\pi}\int_{BZ} \Omega(\mathbf{k})\, d^2k

        an integer (up to the numerical precision set by *nk*) for a group
        of bands that is isolated from the rest of the spectrum by a gap
        everywhere in the Brillouin zone. See *berry_curvature*.

        :param bands: Positive integer, or list of positive integers. Band
            index, or indices of a group of bands (e.g. all occupied bands
            below a gap).
        :param nk: Positive integer, or tuple of 2 positive integers.
            Default value 30. Number of k-points along each reciprocal
            lattice vector.

        :returns:
            * **chern** -- Real number, close to an integer.
        '''
        return self.berry_curvature(bands, nk).sum() / (2*PI)

    def plot_dos(
        self,
        nk: int | tuple[int, int] = 30,
        broadening: float = 0.05,
        kernel: str = 'gaussian',
        e_grid: ArrayLike | None = None,
        fs: float = 20,
        lw: float = 2.,
        figsize: tuple[float, float] | None = None,
    ) -> Figure:
        '''
        Plot the (broadened) density of states, obtained by diagonalizing
        :math:`H(\\mathbf{k})` over a uniform Brillouin-zone mesh -- see
        *tbee.dos.density_of_states*.

        :param nk: Positive integer, or tuple of *dim* positive integers.
            Default value 30. Number of k-points along each reciprocal
            lattice vector.
        :param broadening: Positive real number. Default value 0.05. Kernel width.
        :param kernel: String. Default value 'gaussian'. 'gaussian' or 'lorentzian'.
        :param e_grid: Real ndarray. Default value None. Energies at which to
            evaluate the density of states.
        :param fs: Positive number. Default value 20. Fontsize.
        :param lw: Positive number. Default value 2. Linewidth.
        :param figsize: Tuple. Default value None. Figure size.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.positive_real(fs, 'fs')
        error_handling.positive_real(lw, 'lw')
        error_handling.tuple_2elem(figsize, 'figsize')
        en = self.mesh_bands(nk)
        e_grid, rho = dos.density_of_states(en, e_grid=e_grid,
                                                                   broadening=broadening, kernel=kernel)
        fig, ax = plt.subplots(figsize=figsize)
        ax.plot(e_grid, rho, 'b', lw=lw)
        ax.fill_between(e_grid, rho, color='b', alpha=0.2)
        ax.set_xlim([e_grid[0], e_grid[-1]])
        ax.set_ylim([0., None])
        ax.set_title('Density of states', fontsize=fs)
        ax.set_xlabel('$E$', fontsize=fs)
        ax.set_ylabel(r'$\rho(E)$', fontsize=fs)
        for label in ax.xaxis.get_majorticklabels():
            label.set_fontsize(fs)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def plot_bands(
        self,
        node_labels: list[str] | None = None,
        fs: float = 20,
        lw: float = 2.,
        ms: float = 0.,
        c: str = 'b',
        lims: tuple[float, float] | None = None,
        figsize: tuple[float, float] | None = None,
    ) -> Figure:
        '''
        Plot the band structure computed by *k_path* or *get_bands*.

        :param node_labels: List of strings. Default value None. Labels of the
            high-symmetry points passed to *k_path*.
        :param fs: Positive number. Default value 20. Fontsize.
        :param lw: Positive number. Default value 2. Linewidth.
        :param ms: Positive number. Default value 0. Marker size.
        :param c: Default value 'b'. Line color.
        :param lims: List. Default value None. Energy plot limits.
        :param figsize: Tuple. Default value None. Figure size.

        :returns:
            * **fig** -- Figure.
        '''
        error_handling.empty_ndarray(self.en, 'get_bands or k_path')
        error_handling.positive_real(fs, 'fs')
        error_handling.positive_real(lw, 'lw')
        error_handling.lims(lims)
        error_handling.tuple_2elem(figsize, 'figsize')
        fig, ax = plt.subplots(figsize=figsize)
        for n in range(self.norb):
            ax.plot(self.ks_dist, self.en[:, n], c=c, lw=lw, marker='o', ms=ms)
        for node in self.nodes:
            ax.axvline(node, color='k', lw=0.5)
        ax.set_xlim([self.ks_dist[0], self.ks_dist[-1]])
        if lims is not None:
            ax.set_ylim(lims)
        if node_labels is not None:
            error_handling.ndarray(np.array(node_labels), 'node_labels', len(self.nodes))
            ax.set_xticks(self.nodes)
            ax.set_xticklabels(node_labels, fontsize=fs)
        ax.set_ylabel('$E$', fontsize=fs)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(fs)
        fig.set_layout_engine('tight')
        plt.draw()
        return fig

    def show(self) -> None:
        '''
        Emulate Matplotlib method plt.show().
        '''
        plt.show()


def ribbon(
    lat: Lattice,
    list_hop: list[dict],
    width: int,
    direction: int = 1,
    onsite: dict | None = None,
    spin: bool = False,
) -> KSpace:
    r'''
    Cut a ribbon out of a 2D periodic model: periodic along one primitive
    vector, finite (open boundary, *width* unit cells) along the other.
    This is the standard way to see edge states in a band structure (e.g.
    the zero-energy edge band of a zigzag graphene ribbon, or the helical
    edge states of a Kane-Mele ribbon).

    :param lat: **Lattice** class instance (2D, i.e. two primitive
        vectors). Only *unit_cell* and *prim_vec* are used.
    :param list_hop: List of dictionaries, in the same format passed to
        *KSpace.set_hopping* -- the hoppings of the periodic (2D) model
        that the ribbon is cut from.
    :param width: Positive integer. Number of unit cells across the ribbon.
    :param direction: 0 or 1. Default value 1. Which primitive vector
        (``lat.prim_vec[direction]``) becomes finite; the other stays
        periodic.
    :param onsite: Dictionary. Default value None. Onsite energies, in the
        same format passed to *KSpace.set_onsite* -- applied identically
        on every row of the ribbon.
    :param spin: Boolean. Default value False. See *KSpace*.

    :returns:
        * **rib** -- **KSpace** instance, 1D-periodic, with
          ``width * len(lat.unit_cell)`` sites (each site of *lat*,
          repeated once per row across the ribbon; row *w*'s copy of site
          *i* is orbital ``w*len(lat.unit_cell) + i``).

    Example usage::

        # zigzag graphene ribbon, 20 unit cells wide
        list_hop = [{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                          {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                          {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}]
        rib = ribbon(lat, list_hop, width=20)
    '''
    error_handling.lat(lat)
    error_handling.dim_2(len(lat.prim_vec))
    error_handling.positive_int(width, 'width')
    error_handling.direction(direction)
    periodic = 1 - direction
    n_sites = len(lat.unit_cell)
    a_dir = np.array(lat.prim_vec[direction])
    new_unit_cell = []
    for w in range(width):
        for dic in lat.unit_cell:
            r0 = np.array(dic['r0']) + w*a_dir
            new_unit_cell.append({'tag': dic['tag'], 'r0': (float(r0[0]), float(r0[1]))})
    new_lat = Lattice(unit_cell=new_unit_cell, prim_vec=[lat.prim_vec[periodic]])
    rib = KSpace(new_lat, spin=spin)
    new_list_hop = []
    for dic in list_hop:
        w2_shift = dic['R'][direction]
        n_periodic = dic['R'][periodic]
        for w in range(width):
            w2 = w + w2_shift
            if 0 <= w2 < width:
                new_list_hop.append({'i': w*n_sites + dic['i'],
                                                    'j': w2*n_sites + dic['j'],
                                                    'R': (n_periodic,),
                                                    't': dic['t']})
    rib.set_hopping(new_list_hop)
    if onsite is not None:
        rib.set_onsite(onsite)
    return rib
