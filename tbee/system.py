from __future__ import annotations

from typing import Callable

import numpy as np
from numpy.typing import NDArray
import scipy.sparse as sparse
import scipy.linalg as LA
import numpy.random as rand
import numpy.char as npc
from math import sin, cos
import tbee.error_handling as error_handling
from tbee.lattice import Lattice, COOR_DTYPE


PI = np.pi
ATOL = 1e-3
HOP_DTYPE = [('n', 'u2'), ('i', 'u4'), ('j', 'u4'),
                       ('ang', 'f8'), ('tag', 'U2'), ('t', 'c16')]


class System():
    '''
    Solve the Tight-Binding eigenvalue problem of a lattice defined 
    by the class **lattice**.

    :param lat: **lattice** class instance.
    '''

    def __init__(self, lat: Lattice) -> None:
        error_handling.lat(lat)
        self.lat = lat
        self.sites = self.lat.sites  # used to check if sites changes
        self.coor_hop = np.array([], dtype=COOR_DTYPE)
        self.vec_hop = np.array([], dtype=[('dis', 'f8'),  ('ang', 'f8')]) # Hopping distances and angles
        self.dist_uni = np.array([], 'f8')  # Different hopping distances
        self.store_hop = {}  #  Store the relevant hoppings (dynamic programming)
        self.hop = np.array([], dtype=HOP_DTYPE) #  Hoppings to build-up the Hamiltonian
        self.onsite = np.array([], 'c16')  #  Onsite energies
        self.ham = sparse.csr_matrix(([], ([], [])), shape=(self.lat.sites, self.lat.sites))  # Hamiltonian
        self.en = np.array([], 'c16')  # Eigenenergies
        self.rn = np.array([], 'c16')  # Right eigenvectors: H |rn> = en |rn>
        self.ln = np.array([], 'c16')  # Left eigenvectors:  <ln| H = en <ln|
        self.intensity = np.array([], 'f8')  # Intensities (|rn|**2)
        self.pola = np.array([], 'f8')  # sublattices polarisation (|rn^{(S)}|**2)
        self.petermann = np.array([], 'f8')  # Inverse Participation Ratio
        self.nmax = 0  # number of different hoppings

    def clear_hopping(self) -> None:
        '''
        Clear structured array *hop*.
        '''
        self.hop = np.array([], dtype=HOP_DTYPE)

    def get_distances(self) -> None:
        '''
        Private method.
        Get distances and angles of the edges.
        '''
        error_handling.sites(self.lat.sites)
        dif_x = self.lat.coor['x'] - self.lat.coor['x'].reshape(self.lat.sites, 1)
        dif_y = self.lat.coor['y'] - self.lat.coor['y'].reshape(self.lat.sites, 1)
        dist = np.sqrt(dif_x ** 2 + dif_y ** 2)
        ang = (180 / PI * np.arctan2(dif_y, dif_x))
        self.vec_hop = np.zeros(dist.shape, dtype=[('dis', 'f8'),  ('ang', 'f8')])
        self.vec_hop['dis'] = dist
        self.vec_hop['ang'] = ang
        self.dist_uni = np.unique(self.vec_hop['dis'].round(4))

    def print_distances(self, n: int = 1) -> None:
        r'''
        Print distances and positive angles (in degrees) :math:`\phi_+\in[0, 180)`
        of the nth shortest edges. Negative angles are given by:
        :math:`\phi_-= \phi_+-180` and :math:`\phi_+\in[-180, 0)`.

        :param n: Positive integer. Number of shortest edges.
        '''
        error_handling.sites(self.lat.sites)
        self.get_distances()
        self.nmax = len(self.dist_uni) - 1
        error_handling.positive_int_lim(n, 'n', self.nmax)
        print('\n{} different distances between sites:'.format(self.nmax))
        print('\nDistances between sites:')
        for i, d in enumerate(self.dist_uni[1: n+1]):
            if i == 0:
                hop_name = 'st'
            elif i == 1:
                hop_name = 'nd'
            elif i == 2:
                hop_name = 'rd'
            else:
                hop_name = 'th'
            print('{}{} hopping, length: {:.3f}'.format(i+1, hop_name, d))
            print('\twith positive angles:')
            positive_ang = self.vec_hop['ang'][np.isclose(d, self.vec_hop['dis'], atol=ATOL) &
                                                                   (self.vec_hop['ang'] >= 0.) &
                                                                   (self.vec_hop['ang'] < 180.)]
            print('\t', np.unique(positive_ang.round(4)))

    def set_onsite(self, dict_onsite: dict[str, complex]) -> None:
        '''
        Set onsite energies.

        :param on:  Array. Sublattice onsite energies.

        Example usage::

            # Line-Centered Square lattice
            sys.set_onsite({'a': -1j, 'b': -2j})
        '''
        error_handling.sites(self.lat.sites)
        error_handling.set_onsite(dict_onsite, self.lat.tags)
        self.onsite = np.zeros(self.lat.sites, 'c16')
        for tag, on in dict_onsite.items():
            self.onsite[self.lat.coor['tag'] ==tag] = on

    def fill_store_hop(self, n: int) -> None:
        '''
        Private method.

        Store in *store_hop* indices (with :math:`i < j`), positive angles, and tags
        of a given type of hopping.
        '''
        ind = np.argwhere(np.isclose(self.dist_uni[n], self.vec_hop['dis'], atol=ATOL))
        ind_up = ind[ind[:, 1] > ind[:, 0]]
        hop = np.zeros(len(ind_up), dtype=HOP_DTYPE[:-1])
        hop['i'] = ind_up[:, 0]
        hop['j'] = ind_up[:, 1]
        hop['ang'] = self.vec_hop['ang'][ind_up[:, 0], ind_up[:, 1]]
        hop['tag'] = npc.add(self.lat.coor['tag'][ind_up[:, 0]], 
                                         self.lat.coor['tag'][ind_up[:, 1]])
        self.store_hop[n] = hop

    def set_hopping(self, list_hop: list[dict], upper_part: bool = True) -> None:
        r'''
        Set lattice hoppings.

        :param list_hop: List of Dictionaries.
            Dictionary with keys ('n', 'ang', 'tag', 't') where:

                * 'n' Positive integer, type of hoppings:

                    * 'n': 1 for nearest neighbours.
                    * 'n': 2 for next-nearest neighbours.  
                    * 'n': 3 for next-next-nearest neighbours.  
                    * etc...

                * 'ang' value, float, angle, in deg, of the hoppings. (optional).

                    Hopping angles are given by the method *print_distances*.

                        * If :math:`ang \in[0, 180)`, fill the Hamiltonian upper part.
                        * If :math:`ang \in[-180, 0)`, fill the Hamiltonian lower part.

                * 'tag' string of length 2  (optional).

                    Hopping tags.

                * 't' Complex number.

                    Hopping value.

        :param upper_part: Boolean. Default value True.

            * True get hoppings with (:math:`i<j`) *i.e.* fill the Hamiltonian upper part.
            * False get hoppings with (:math:`i>j`) *i.e.* fill the Hamiltonian lower part.

        Example usage::

            # fill upper part:
            sys.set_hopping([{'n': 1, t: 1.}])
            # fill lower part:
            sys.set_hopping([{'n': 1, t: 1.}], upper_part=False)
            # fill upper part: specifying the angles:
            sys.set_hopping([{'n': 1, 'ang': 0., t: 1.}, {'n': 1, 'ang': 90,  t: 2.}])
            # fill lower part:
            sys.set_hopping([{'n': 1, 'ang': -180., t: 1.}, {'n': 1, 'ang': -90,  t: 2.}], upper_part=False)
            # fill upper part: specifying the tags:
            sys.set_hopping([{'n': 1, 'tag': 'ab', t: 1.}, {'n': 1, 'tag': 'ba',  t: 2.}])
            # fill lower part:
            sys.set_hopping([{'n': 1, 'tag': 'ab', t: 1.}, {'n': 1, 'tag': 'ba',  t: 2.}], upper_part=False)
            # fill upper part: specifying the angles and tags:
            sys.set_hopping([{'n': 1, 'ang': 0., 'tag': 'ab', t: 1.}, 
                                        {'n': 1, 'ang': 0., 'tag': 'ba',  t: 2.},
                                        {'n': 1, 'ang': 90., 'tag': 'ab', t: 3.}, 
                                        {'n': 1, 'ang': 90., 'tag': 'ba',  t: 4.}])
            # fill lower part:
            sys.set_hopping([{'n': 1, 'ang': 0., 'tag': 'ab', t: 1.}, 
                                        {'n': 1, 'ang': 0., 'tag': 'ba',  t: 2.},
                                        {'n': 1, 'ang': 90., 'tag': 'ab', t: 3.}, 
                                        {'n': 1, 'ang': 90., 'tag': 'ba',  t: 4.}]), upper_part=False)

        .. note::

            A Hermitian hopping matrix can be build-up only using 
            its upper part OR only using its lower part. The full matrix is then
            automatic built by Hermitian conjugaison. 
            
            If both upper AND lower parts are used to build up the hopping matrix.
            non Hermitian conjugaison is not performed *i.e.* non-Hermitian hopping matrix
            can be built.
        '''
        error_handling.sites(self.lat.sites)
        error_handling.boolean(upper_part, 'upper_part')
        self.get_distances()
        self.nmax = len(self.dist_uni) - 1
        error_handling.set_hopping(list_hop, self.nmax)
        list_n = np.unique([dic['n'] for dic in list_hop])
        # fill, if needed self.store_hop
        self.check_sites()
        for n in list_n:
            if n not in self.store_hop:
                self.fill_store_hop(n)
        # fill self.hop
        for dic in list_hop:
            if len(dic) == 2:
                size = len(self.store_hop[dic['n']])
                if upper_part:
                    mask = (self.hop['n'] == dic['n']) & (self.hop['i'] < self.hop['j'])
                else:
                    mask = (self.hop['n'] == dic['n']) & (self.hop['i'] > self.hop['j'])
                if np.sum(mask):
                    self.hop = self.hop[np.logical_not(mask)]
                ind = np.ones(size, bool)
                hop = self.set_given_hopping(dic['n'], size, dic, ind, upper_part=upper_part)
            elif len(dic) == 3 and 'ang' in dic:
                error_handling.angle(dic['ang'], np.unique(self.store_hop[dic['n']]['ang']), upper_part)
                if dic['ang'] >= 0:
                    ang_store = dic['ang']
                else:
                    ang_store = dic['ang'] + 180.
                size = np.sum(np.isclose(ang_store, self.store_hop[dic['n']]['ang'], atol=ATOL))
                mask = (self.hop['n'] == dic['n']) & np.isclose(self.hop['ang'], dic['ang'], atol=ATOL)
                if np.sum(mask):
                    self.hop = self.hop[np.logical_not(mask)]
                ind = np.isclose(ang_store, self.store_hop[dic['n']]['ang'], atol=ATOL)
                error_handling.index(ind, dic)
                hop = self.set_given_hopping(dic['n'], size, dic, ind, upper_part=upper_part)
            elif len(dic) == 3 and 'tag' in dic:
                if upper_part:
                    tag_store = dic['tag']
                else:
                    tag_store = dic['tag'][::-1]
                size = np.sum(self.store_hop[dic['n']]['tag'] == tag_store)
                mask = (self.hop['n'] == dic['n']) & (self.hop['tag'] == dic['tag'])
                if upper_part:
                    mask = self.hop['n'] == dic['n'] & (self.hop['tag'] == dic['tag']) & (self.hop['i'] < self.hop['j'])
                else:
                    mask = self.hop['n'] == dic['n'] & (self.hop['tag'] == dic['tag']) & (self.hop['i'] > self.hop['j'])
                if np.sum(mask):
                    self.hop = self.hop[np.logical_not(mask)]
                ind = self.store_hop[dic['n']]['tag'] == tag_store
                error_handling.index(ind, dic)
                hop = self.set_given_hopping(dic['n'], size, dic, ind, upper_part=upper_part)
            else:
                error_handling.angle(dic['ang'], np.unique(self.store_hop[dic['n']]['ang']), upper_part=upper_part)
                error_handling.tag(dic['tag'], np.unique(self.store_hop[dic['n']]['tag']))
                if dic['ang'] >= 0:
                    ang_store = dic['ang']
                else:
                    ang_store = dic['ang'] + 180.
                if upper_part:
                    tag_store = dic['tag']
                else:
                    tag_store = dic['tag'][::-1]
                size = np.sum((self.store_hop[dic['n']]['tag'] == tag_store) & 
                                       (np.isclose(ang_store, self.store_hop[dic['n']]['ang'], atol=ATOL)))
                bool1 = (self.hop['n'] == dic['n']) & (self.hop['tag'] == dic['tag'])
                bool2 = np.isclose(self.hop['ang'], dic['ang'], atol=ATOL)
                mask = bool1 & bool2
                if np.sum(mask):
                    self.hop = self.hop[np.logical_not(mask)]
                ind = ((self.store_hop[dic['n']]['tag'] == tag_store) & 
                          (np.isclose(ang_store, self.store_hop[dic['n']]['ang'], atol=1)))
                error_handling.index(ind, dic)
                hop = self.set_given_hopping(dic['n'], size, dic, ind, upper_part=upper_part)
            self.hop = np.concatenate([self.hop, hop])

    def check_sites(self) -> None:
        '''
        Private method.
        Check if the number of sites was changed after calling the 
        method system.set_hopping().
        '''
        if self.sites != self.lat.sites:
             self.store_hop = {}
             self.sites = self.lat.sites

    def set_given_hopping(
        self, n: int, size: int, dic: dict, mask: NDArray, upper_part: bool,
    ) -> NDArray:
        '''
        Private method.
        Fill self.hop. 

        :param n: Integer. Hopping type.
        :param size: Integer. Number of hoppings.
        :param doc: Dictionary. Hopping dictionary.
        :param mask: np.ndarray. Mask.
        :param upper_part: Boolean. If True, self.hop['i'] < self.hop['j'].
        '''
        hop = np.empty(size, dtype=HOP_DTYPE)
        hop['n'] = dic['n']
        hop['t'] = dic['t']
        if upper_part:
            hop['i'] = self.store_hop[n]['i'][mask]
            hop['j'] = self.store_hop[n]['j'][mask]
            hop['ang'] = self.store_hop[n]['ang'][mask]
            hop['tag'] = self.store_hop[n]['tag'][mask]
        else:
            hop['i'] = self.store_hop[n]['j'][mask]
            hop['j'] = self.store_hop[n]['i'][mask]
            hop['ang'] = self.store_hop[n]['ang'][mask] - 180
            hop['tag'] = npc.add(self.lat.coor['tag'][hop['i']], 
                                            self.lat.coor['tag'][hop['j']])
        return hop

    def set_hopping_manual(self, dict_hop: dict[tuple[int, int], complex], upper_part: bool = True) -> None:
        '''
        Set hoppings manually.

        :param dict_hop: Dictionary of hoppings.
            key: hopping indices, val: hopping values.

        :parameter upper_part: Boolean. 

            * True, fill the Hamiltonian upper part.
            * False, fill the Hamiltonian lower part.  
        '''
        hop = np.zeros(len(dict_hop), dtype=HOP_DTYPE)
        i = [h[0] for h in dict_hop.keys()]
        j = [h[1] for h in dict_hop.keys()]
        t = [val for val in dict_hop.values()]
        hop['i'],  hop['j']= i, j
        hop['t'] = t 
        hop['tag'] = npc.add(self.lat.coor['tag'][i], 
                                        self.lat.coor['tag'][j])
        ang = 180 / PI * np.arctan2(self.lat.coor['y'][j]-self.lat.coor['y'][i],
                                                    self.lat.coor['x'][j]-self.lat.coor['x'][i])
        if upper_part:
            ang[ang < 0] += 180
        else:
            ang[ang >= 0] -= 180
        hop['ang'] = ang
        self.hop = np.concatenate([self.hop, hop])

    def set_hopping_dis(self, alpha: complex) -> None:
        '''
        Set uniform hopping disorder. 

        :param alpha: Complex or Real number. Disorder stength.

        Example usage::

            sys.set_hopping_dis(alpha=0.1)

        '''
        error_handling.empty_hop(self.hop)
        error_handling.number(alpha, 'alpha')
        self.hop['t'] *= 1. + alpha * rand.uniform(-1., 1., len(self.hop))

    def set_peierls_phase(self, phase: Callable[..., NDArray]) -> None:
        r'''
        Apply the Peierls substitution to the existing hoppings, to
        capture the effect of an orbital magnetic field:

        .. math::

            t_{ij} \to t_{ij}\, e^{i\phi_{ij}}\, ,\quad
            \phi_{ij} = \frac{2\pi}{\Phi_0}\int_{\mathbf{r}_i}^{\mathbf{r}_j}
            \mathbf{A}\cdot d\mathbf{l}

        where :math:`\mathbf{A}` is the vector potential, integrated along
        the straight bond from site :math:`i` to site :math:`j`.

        Must be called after *set_hopping* / *set_hopping_manual* (it
        rescales the existing hoppings *in place*) and before *get_ham*.
        For a uniform perpendicular field, use the convenience method
        *set_magnetic_field* instead.

        :param phase: Callable. ``phase(xi, yi, xj, yj)`` returns
            :math:`\phi_{ij}`, the (real-valued) Peierls phase for the bond
            from :math:`(x_i, y_i)` to :math:`(x_j, y_j)`. Called with
            Numpy arrays (one value per hopping in *sys.hop*).

        .. note::

            *get_ham* automatically assigns the reversed bond its complex
            conjugate, so the Hamiltonian stays Hermitian as long as
            *phase* is antisymmetric under swapping :math:`i` and
            :math:`j` -- true for the line integral of any vector
            potential, since reversing the integration path negates it.

        Example usage::

            # Peierls phase from a uniform field via the Landau gauge
            # A = (0, B x): captures the same physics as set_magnetic_field,
            # just in a different (equally valid) gauge.
            B = 0.05
            sys.set_peierls_phase(lambda xi, yi, xj, yj: B * (xj - xi) * (xi + xj) / 2)
        '''
        error_handling.empty_hop(self.hop)
        error_handling.is_callable(phase, 'phase')
        xi = self.lat.coor['x'][self.hop['i']]
        yi = self.lat.coor['y'][self.hop['i']]
        xj = self.lat.coor['x'][self.hop['j']]
        yj = self.lat.coor['y'][self.hop['j']]
        self.hop['t'] = self.hop['t'] * np.exp(1j * phase(xi, yi, xj, yj))

    def set_magnetic_field(self, alpha: float) -> None:
        r'''
        Set a uniform perpendicular magnetic field via the Peierls
        substitution (see *set_peierls_phase*), using the symmetric gauge
        :math:`\mathbf{A} = \frac{B}{2}(-y, x)`:

        .. math::

            \phi_{ij} = \pi\alpha\,(x_iy_j - x_jy_i)

        :param alpha: Real number. Flux density :math:`B/\Phi_0`, in flux
            quanta per unit area (in the lattice's length units) -- *i.e.*
            the flux through a region of area :math:`S` is
            :math:`\alpha S` flux quanta.

        Example usage::

            # one flux quantum per 100 unit cells of a lattice with
            # lattice constant 1:
            sys.set_magnetic_field(alpha=0.01)
        '''
        error_handling.real_number(alpha, 'alpha')
        self.set_peierls_phase(lambda xi, yi, xj, yj: PI * alpha * (xi*yj - xj*yi))

    def set_onsite_dis(self, alpha: complex) -> None:
        '''
        Set uniform onsite disorder. 

        :param alpha: Complex or Real number. Disorder stength.

        Example usage::

            sys.set_onsite_dis(alpha=0.1)

        '''
        error_handling.empty_onsite(self.onsite)
        error_handling.number(alpha, 'alpha')
        self.onsite += alpha * rand.uniform(-1., 1., self.lat.sites)

    def set_onsite_def(self, onsite_def: dict[int, complex]) -> None:
        '''
        Set specific onsite energies.

        :param onsite_def:  Dictionary. 
            key: site indices, val: onsite values.

        Example usage::

            set_onsite_def(0: 1., 1: -1j)
        '''
        error_handling.empty_onsite(self.onsite)
        error_handling.set_onsite_def(onsite_def, self.lat.sites)
        for i, o in onsite_def.items():
            self.onsite[i] = o

    def set_hopping_def(self, hopping_def: dict[tuple[int, int], complex]) -> None:
        '''
        Set specific hoppings. 

        :param hopping_def:  Dictionary of hoppings. 
            key: hopping indices, val: hopping values. 

        Example usage::

            sys.set_hopping_def({(0, 1): 1., (1, 2): -1j})
        '''
        error_handling.empty_hop(self.hop)
        error_handling.set_hopping_def(self.hop, hopping_def, self.lat.sites)
        for key, val in hopping_def.items():
            cond = (self.hop['i'] == key[0]) & (self.hop['j'] == key[1])
            self.hop['t'][cond] = val
            self.hop['ang'] = self.vec_hop['ang'][key[0], key[1]]
            self.hop['tag'] = npc.add(self.lat.coor['tag'][key[0]],
                                                   self.lat.coor['tag'][key[1]])

    def set_new_hopping(self, list_hop: list[dict], ind: NDArray) -> None:
        '''
        Private method.
        Set new hoppings.

        :param list_hop: List of Dictionary (see set_hopping definition).
        :param ind: List. List of indices.
        '''
        for dic in list_hop:
            if len(dic) == 2:
                self.hop['t'][ind] = dic['t']
            elif len(dic) == 3 and 'ang' in dic:
                self.hop['t'][ind & (self.hop['ang'] == dic['ang'])] = dic['t']
            elif len(dic) == 3 and 'tag' in dic:
                self.hop['t'][ind & (self.hop['tag'] == dic['tag'])] = dic['t']
            else:
                self.hop['t'][ind & (self.hop['tag'] == dic['tag'])
                                        & (self.hop['ang'] == dic['ang'])] = dic['t']

    def find_square(self, xlims: tuple[float, float], ylims: tuple[float, float]) -> NDArray:
        '''
        Private method.
        Find hoppings within the square.

        :param xlims: List or Tuple. :math:`x` interval.
        :param ylims: List or Tuple. :math:`y` interval.
        '''
        error_handling.lims(xlims)
        error_handling.lims(ylims)
        in1 = (self.lat.coor['x'][self.hop['i']] >= xlims[0]) & \
                 (self.lat.coor['y'][self.hop['i']] >= ylims[0]) & \
                 (self.lat.coor['x'][self.hop['j']] >= xlims[0]) & \
                 (self.lat.coor['y'][self.hop['j']] >= ylims[0])
        in2 = (self.lat.coor['x'][self.hop['i']] <= xlims[1]) & \
                 (self.lat.coor['y'][self.hop['i']] <= ylims[1]) & \
                 (self.lat.coor['x'][self.hop['j']] <= xlims[1]) & \
                 (self.lat.coor['y'][self.hop['j']] <= ylims[1])
        return in1 * in2

    def find_ellipse(self, rx: float, ry: float, x0: float, y0: float) -> NDArray:
        '''
        Private method.
        Find hoppings within the ellipse.

        :param rx: Positive Float. Radius along :math:`x`. 
        :param ry: Positive Float. Radius along :math:`y`.
        :param x0: Float. Defalut value 0. :math:`x` center. 
        :param y0: Float. Defalut value 0. :math:`x` center.
        '''
        in1 = (self.lat.coor['x'][self.hop['i']] - x0) ** 2 / rx ** 2 + \
                 (self.lat.coor['y'][self.hop['i']] - y0) ** 2 / ry ** 2 <= 1.
        in2 = (self.lat.coor['x'][self.hop['j']] - x0) ** 2 / rx ** 2 + \
                 (self.lat.coor['y'][self.hop['j']] - y0) ** 2 / ry ** 2 <= 1.
        return in1 * in2

    def change_hopping_square(
        self, list_hop: list[dict], xlims: tuple[float, float], ylims: tuple[float, float] = [-1., 1.],
    ) -> None:
        '''
        Change hopping values.

        :param list_hop: List of Dictionary (see set_hopping definition).
        :param xlims: List or Tuple. :math:`x` interval.
        :param ylims: List or Tuple. :math:`y` interval.
        '''
        error_handling.empty_hop(self.hop)
        error_handling.set_hopping(list_hop, self.nmax)
        ind = self.find_square(xlims, ylims)
        self.set_new_hopping(list_hop, ind)

    def change_hopping_ellipse(
        self, list_hop: list[dict], rx: float, ry: float, x0: float = 0., y0: float = 0.,
    ) -> None:
        '''
        Change hopping values.

        :param list_hop: List of Dictionary (see set_hopping definition).
        :param rx: Positive Float. Radius along :math:`x`. 
        :param ry: Positive Float. Radius along :math:`y`.
        :param x0: Float. Default value 0. :math:`x` center. 
        :param y0: Float. Default value 0. :math:`y` center.
        '''
        error_handling.empty_hop(self.hop)
        error_handling.set_hopping(list_hop, self.nmax)
        error_handling.positive_real(rx, 'rx')
        error_handling.positive_real(ry, 'rx')
        error_handling.real_number(x0, 'x0')
        error_handling.real_number(y0, 'y0')
        ind = self.find_ellipse(rx, ry, x0, y0)
        self.set_new_hopping(list_hop, ind)

    def get_coor_hop(self) -> None:
        '''
        Get the site coordinates in hopping space
        only considering the nearest neighbours hoppings.
        '''
        error_handling.empty_hop(self.hop)
        visited = np.zeros(self.lat.sites, 'u2')
        self.coor_hop = np.zeros(self.lat.sites, dtype=COOR_DTYPE)
        self.coor_hop['tag'] = self.lat.coor['tag']
        hop = self.hop[self.hop['n'] == 1]
        hop_down = np.copy(hop)
        hop_down['i'] = hop['j']
        hop_down['j'] = hop[ 'i']
        hop_down['ang'] = -180 + hop['ang']
        hop = np.concatenate([hop, hop_down])
        i_visit = np.min(hop['i'])
        while True:
            hs = hop[hop['i'] == i_visit]
            for h in hs:
                if visited[h['j']] == 2:
                    continue
                self.coor_hop['x'][h['j']] = self.coor_hop['x'][i_visit] + \
                    h['t'].real*cos(PI / 180 * h['ang'])
                self.coor_hop['y'][h['j']] = self.coor_hop['y'][i_visit] + \
                    h['t'].real*sin(PI / 180 * h['ang'])
                visited[h['j']] = 1
            visited[i_visit] = 2
            explored = np.argwhere(visited == 1)
            if not explored.any():
                break
            i_visit = explored[0, 0]

    def get_ham(self) -> None:
        '''
        Get the Tight-Binding Hamiltonian using sys.hop.
        '''
        error_handling.empty_hop(self.hop)
        error_handling.hop_sites(self.hop, self.lat.sites)
        if np.all(self.hop['ang'] >= 0) or np.all(self.hop['ang'] < 0):
            self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                            shape=(self.lat.sites, self.lat.sites)) \
                           + sparse.csr_matrix((self.hop['t'].conj(), (self.hop['j'], self.hop['i'])), 
                                                            shape=(self.lat.sites, self.lat.sites))
        else:
            self.ham = sparse.csr_matrix((self.hop['t'], (self.hop['i'], self.hop['j'])), 
                                                            shape=(self.lat.sites, self.lat.sites))
        if self.onsite.size == self.lat.sites:
            self.ham += sparse.diags(self.onsite, 0)

    def get_eig(self, eigenvec: bool = False, left: bool = False) -> None:
        '''
        Get the eigenergies, eigenvectors and polarisation.

        :param eigenvec: Boolean. Default value False. 
            If True, get the eigenvectors.
        :param left: Boolean. Default value False. 
            If True, get the left eigenvectors too. 
            Relevant for non-Hermitian matrices.
        '''
        error_handling.empty_ham(self.ham)
        error_handling.boolean(eigenvec, 'eigenvec')
        error_handling.boolean(left, 'left')
        if eigenvec:
            if (self.ham.conj().T != self.ham).nnz:
                if not left:
                    self.en, self.rn = LA.eig(self.ham.toarray())
                else:
                    self.en, self.rn, self.ln = LA.eig(self.ham.toarray(), left=left)
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
                self.rn = self.rn[:, ind]
                if self.ln.size:
                    self.ln = self.ln[:, ind]
            else:
                self.en, self.rn = LA.eigh(self.ham.toarray())
            self.intensity = np.abs(self.rn) ** 2
            self.pola = np.zeros((self.lat.sites, len(self.lat.tags)))
            for i, tag in enumerate(self.lat.tags):
                self.pola[:, i] = np.sum(self.intensity[self.lat.coor['tag'] == tag, :], axis=0)
        else:
            if (self.ham.conj().T != self.ham).nnz:
                self.en = LA.eigvals(self.ham.toarray())
                ind = np.argsort(self.en.real)
                self.en = self.en[ind]
            else:
                self.en = LA.eigvalsh(self.ham.toarray())

    def get_ipr(self) -> None:
        r'''
        Get the Inverse Participation Ratio: 

        .. math:: 

            IPR_n = |\sum_i\psi_i^{n}|^4\, .
        '''
        error_handling.empty_ndarray(self.rn, 'sys.get_eig(eigenvec=True)')
        self.ipr = np.sum(self.intensity ** 2, axis=0)

    def get_petermann(self) -> None:
        r'''
        Get the Petermann factor: 
        
        .. math::

            K_n = \frac{\langle\psi_L^{n}|\psi_L^{n}\rangle\langle\psi_R^{n}|\psi_R^{n}\rangle}{\langle\psi_L^{n}|\psi_R^{n}\rangle}\, .

        .. note::

            LA.eig fixes the norm such that :math:`\langle\psi_L^{n}|\psi_L^{n}\rangle = 1` and :math:`\langle\psi_R^{n}|\psi_R^{n}\rangle = 1`.
        '''
        if not (self.ham.conj().T != self.ham).nnz:
            self.petermann = np.ones(self.lat.sites)
            return
        error_handling.empty_ndarray(self.ln, 'sys.get_eig(eigenvec=True, left=True)')
        left_right = np.sum(self.ln * np.conjugate(self.rn), axis=0).real
        self.petermann = 1. / left_right ** 2

    def get_intensity_pola_max(self, tag_pola: str) -> NDArray[np.float64]:
        '''
        Get the state with largest polarization on one sublattice.

        :param tag_pola: One-character string. Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        error_handling.empty_ndarray(self.rn, 'sys.get_eig(eigenvec=True)')
        error_handling.tag(tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmax(self.pola[:, i_tag])
        print('State with polarization: {:.5f}'.format(self.pola[ind, i_tag].item()))
        return self.intensity[:, ind]

    def get_intensity_pola_min(self, tag_pola: str) -> NDArray[np.float64]:
        '''
        Get the state with smallest polarization on one sublattice.

        :param tag_pola: One-character string. Sublattice tag.

        :returns:
            * **intensity** -- Intensity of max polarized state on *tag*.
        '''
        error_handling.empty_ndarray(self.rn, 'sys.get_eig(eigenvec=True)')
        error_handling.tag(tag_pola, self.lat.tags)
        i_tag = self.lat.tags == tag_pola
        ind = np.argmin(self.pola[:, i_tag])
        print('State with polarization: {:.5f}'.format(self.pola[ind, i_tag].item()))
        return self.intensity[:, ind]

    def get_intensity_en(self, lims: tuple[float, float]) -> NDArray[np.float64]:
        '''
        Get, if any, the intensity of the sum of the states 
        between *lims[0]* and *lims[1]*.

        :param lims: List. lims[0] energy min, lims[1] energy max.

        :returns:
            * **intensity** -- Sum of the intensities between (lims[0], lims[1]).
        '''
        error_handling.empty_ndarray(self.rn, 'sys.get_eig(eigenvec=True)')
        error_handling.lims(lims)
        ind = np.where((self.en > lims[0]) & (self.en < lims[1]))
        ind = np.ravel(ind)
        print('{} states between {} and {}'.format(len(ind), lims[0], lims[1]))
        return np.sum(self.intensity[:, ind], axis=1)


# Backward-compatible lowercase alias (pre-0.2 API).
system = System
