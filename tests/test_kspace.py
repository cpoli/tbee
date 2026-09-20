from tbee.lattice import lattice
from tbee.kspace import KSpace, reciprocal_vectors, PAULI, ribbon
from tbee.system import system
import unittest
import numpy as np
from math import sqrt


DX, DY = 0.5 * sqrt(3), 0.5


def graphene():
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
    prim_vec = [(2*DX, 0.), (DX, 1.5)]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    kag = KSpace(lat)
    kag.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}])
    return kag, prim_vec


def haldane(t1=1., t2=0.2, phi=np.pi/2, M=0.):
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
    prim_vec = [(2*DX, 0.), (DX, 1.5)]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    hal = KSpace(lat)
    hal.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': t1}])
    for R in [(1, 0), (0, 1), (1, -1)]:
        hal.set_hopping([{'i': 0, 'j': 0, 'R': R, 't': t2*np.exp(1j*phi)}])
        hal.set_hopping([{'i': 1, 'j': 1, 'R': R, 't': t2*np.exp(-1j*phi)}])
    hal.set_onsite({'a': M, 'b': -M})
    return hal


class TestKSpace(unittest.TestCase):
    '''
    Unittest of class **KSpace**.
    '''
    def test_init(self):
        self.assertRaises(TypeError, KSpace, lat=0)

    def test_reciprocal_vectors_2d(self):
        prim_vec = [(2*DX, 0.), (DX, 1.5)]
        a1, a2 = (np.array(v) for v in prim_vec)
        b1, b2 = (np.array(v) for v in reciprocal_vectors(prim_vec))
        self.assertTrue(np.isclose(a1 @ b1, 2*np.pi))
        self.assertTrue(np.isclose(a1 @ b2, 0.))
        self.assertTrue(np.isclose(a2 @ b1, 0.))
        self.assertTrue(np.isclose(a2 @ b2, 2*np.pi))

    def test_reciprocal_vectors_1d(self):
        prim_vec = [(2., 0.)]
        a1 = np.array(prim_vec[0])
        b1, = (np.array(v) for v in reciprocal_vectors(prim_vec))
        self.assertTrue(np.isclose(a1 @ b1, 2*np.pi))

    def test_set_onsite(self):
        kag, _ = graphene()
        self.assertRaises(TypeError, kag.set_onsite, 0)
        self.assertRaises(ValueError, kag.set_onsite, {'z': 1.})
        kag.set_onsite({'a': 1., 'b': -1.})
        self.assertTrue(np.allclose(kag.onsite, [1., -1.]))

    def test_set_hopping_errors(self):
        kag, _ = graphene()
        self.assertRaises(TypeError, kag.set_hopping, 0)
        self.assertRaises(KeyError, kag.set_hopping, [{'i': 0, 'j': 1}])
        self.assertRaises(ValueError, kag.set_hopping, [{'i': 0, 'j': 5, 'R': (0, 0), 't': 1.}])
        self.assertRaises(ValueError, kag.set_hopping, [{'i': 0, 'j': 0, 'R': (0, 0), 't': 1.}])
        self.assertRaises(ValueError, kag.set_hopping, [{'i': 0, 'j': 1, 'R': (0,), 't': 1.}])

    def test_hermitian(self):
        kag, _ = graphene()
        for k in ([0., 0.], [0.3, -1.2], [2.1, 0.7]):
            ham = kag.get_ham(k)
            self.assertTrue(np.allclose(ham, ham.conj().T))

    def test_graphene_dirac_point(self):
        kag, prim_vec = graphene()
        b1, b2 = (np.array(v) for v in reciprocal_vectors(prim_vec))
        K = (b1 - b2) / 3
        en = np.linalg.eigvalsh(kag.get_ham(K))
        self.assertTrue(np.allclose(en, [0., 0.], atol=1e-8))

    def test_graphene_bandwidth(self):
        kag, _ = graphene()
        ks = [[kx, ky] for kx in np.linspace(-4, 4, 15) for ky in np.linspace(-4, 4, 15)]
        en = kag.get_bands(ks)
        self.assertTrue(np.isclose(en.max(), 3., atol=1e-6))
        self.assertTrue(np.isclose(en.min(), -3., atol=1e-6))

    def test_k_path_nodes(self):
        kag, _ = graphene()
        points = [[0., 0.], [1., 0.], [1., 1.]]
        ks_dist, en = kag.k_path(points, nk=5)
        self.assertEqual(len(ks_dist), 2*5 + 1)
        self.assertEqual(en.shape, (len(ks_dist), 2))
        self.assertTrue(np.isclose(kag.nodes[0], 0.))
        self.assertTrue(np.isclose(kag.nodes[-1], ks_dist[-1]))

    def test_clear_hopping(self):
        kag, _ = graphene()
        kag.clear_hopping()
        ham = kag.get_ham([0., 0.])
        self.assertTrue(np.allclose(ham, 0.))

    def test_get_bands_eigenvec(self):
        kag, _ = graphene()
        ks = [[0., 0.], [0.3, -1.2]]
        en, vn = kag.get_bands(ks, eigenvec=True)
        self.assertEqual(en.shape, (2, 2))
        self.assertEqual(vn.shape, (2, 2, 2))
        for i, k in enumerate(ks):
            ham = kag.get_ham(k)
            for n in range(2):
                self.assertTrue(np.allclose(ham @ vn[i, :, n], en[i, n] * vn[i, :, n], atol=1e-8))

    def test_plot_bands_requires_bands(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}],
                              prim_vec=[(2*DX, 0.), (DX, 1.5)])
        kag = KSpace(lat)
        self.assertRaises(RuntimeError, kag.plot_bands)

    def test_plot_bands(self):
        kag, prim_vec = graphene()
        b1, b2 = (np.array(v) for v in reciprocal_vectors(prim_vec))
        Gamma, K, M = np.zeros(2), (b1 - b2) / 3, b1 / 2
        kag.k_path([Gamma, K, M, Gamma], nk=10)
        fig = kag.plot_bands()
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = kag.plot_bands(node_labels=['G', 'K', 'M', 'G'], lims=[-4., 4.])
        self.assertEqual(fig2.__class__.__name__, 'Figure')
        kag.show()

    def test_1d_chain_bands(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat)
        chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
        en = chain.get_ham((0.,))
        self.assertTrue(np.isclose(en[0, 0].real, 2.))
        en_pi = chain.get_ham((np.pi,))
        self.assertTrue(np.isclose(en_pi[0, 0].real, -2.))

    def test_mesh_bands_2d(self):
        kag, _ = graphene()
        en = kag.mesh_bands(nk=10)
        self.assertEqual(en.shape, (100, 2))
        en2 = kag.mesh_bands(nk=(5, 7))
        self.assertEqual(en2.shape, (35, 2))
        self.assertRaises(TypeError, kag.mesh_bands, 1.5)
        self.assertRaises(ValueError, kag.mesh_bands, (5,))

    def test_mesh_bands_1d(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat)
        chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
        en = chain.mesh_bands(nk=40)
        self.assertEqual(en.shape, (40, 1))
        self.assertTrue(np.isclose(en.max(), 2.))
        self.assertTrue(np.isclose(en.min(), -2., atol=1e-1))

    def test_plot_dos(self):
        kag, _ = graphene()
        fig = kag.plot_dos(nk=20)
        self.assertEqual(fig.__class__.__name__, 'Figure')
        fig2 = kag.plot_dos(nk=20, kernel='lorentzian', broadening=0.1)
        self.assertEqual(fig2.__class__.__name__, 'Figure')


class TestTopology(unittest.TestCase):
    '''
    Unittest of Berry curvature / Chern number (class **KSpace**).
    '''
    def test_errors(self):
        kag, _ = graphene()
        self.assertRaises(TypeError, kag.berry_curvature, 1.5)
        self.assertRaises(TypeError, kag.berry_curvature, [0.5])
        self.assertRaises(ValueError, kag.berry_curvature, [0, 0])
        self.assertRaises(ValueError, kag.berry_curvature, [5])
        self.assertRaises(TypeError, kag.berry_curvature, [0], nk=1.5)
        lat_1d = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat_1d)
        chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
        self.assertRaises(ValueError, chain.berry_curvature, [0])

    def test_haldane_topological_phase(self):
        hal = haldane(t1=1., t2=0.2, phi=np.pi/2, M=0.)
        en = hal.mesh_bands(nk=30)
        self.assertTrue(en[:, 1].min() > en[:, 0].max())  # gapped
        chern = hal.chern_number(bands=[0], nk=40)
        self.assertTrue(np.isclose(chern, 1., atol=1e-2))
        # a bare band index (not wrapped in a list) must give the same result.
        self.assertEqual(hal.chern_number(bands=0, nk=40), chern)
        # Berry curvature must sum (not just each plaquette individually)
        # to 2*pi*chern.
        curv = hal.berry_curvature(bands=[0], nk=40)
        self.assertTrue(np.isclose(curv.sum() / (2*np.pi), chern))
        # the two bands carry opposite Chern number (total must vanish,
        # since the full Hamiltonian is a smooth function of k on a closed
        # manifold).
        chern_both = hal.chern_number(bands=[0, 1], nk=40)
        self.assertTrue(np.isclose(chern_both, 0., atol=1e-6))

    def test_haldane_trivial_phase(self):
        # mass term M dominates the topological term -> trivial insulator.
        hal = haldane(t1=1., t2=0.2, phi=np.pi/2, M=2.)
        chern = hal.chern_number(bands=[0], nk=40)
        self.assertTrue(np.isclose(chern, 0., atol=1e-2))

    def test_time_reversal_symmetric_is_trivial(self):
        # t2=0: ordinary (gapped, via M) graphene, time-reversal symmetric
        # -> must be topologically trivial.
        hal = haldane(t1=1., t2=0., phi=np.pi/2, M=0.1)
        chern = hal.chern_number(bands=[0], nk=40)
        self.assertTrue(np.isclose(chern, 0., atol=1e-6))


class TestSpin(unittest.TestCase):
    '''
    Unittest of the spinful mode (``KSpace(lat, spin=True)``).
    '''
    def test_norb_doubles(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat, spin=True)
        self.assertEqual(chain.n_sites, 1)
        self.assertEqual(chain.norb, 2)

    def test_set_hopping_errors(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat, spin=True)
        self.assertRaises(TypeError, chain.set_hopping, [{'i': 0, 'j': 0, 'R': (1,), 't': [1., 2.]}])
        self.assertRaises(TypeError, chain.set_hopping,
                                    [{'i': 0, 'j': 0, 'R': (1,), 't': np.eye(3)}])

    def test_set_onsite_errors(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat, spin=True)
        self.assertRaises(TypeError, chain.set_onsite, {'a': 'z'})
        self.assertRaises(TypeError, chain.set_onsite, {'a': (1., 'z')})
        self.assertRaises(TypeError, chain.set_onsite, {'a': (1., 2., 3.)})

    def test_spin_independent_hopping_is_doubly_degenerate(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat, spin=True)
        chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
        en = np.linalg.eigvalsh(chain.get_ham((0.3,)))
        self.assertTrue(np.isclose(en[0], en[1]))

    def test_onsite_zeeman_splitting(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        chain = KSpace(lat, spin=True)
        chain.set_hopping([{'i': 0, 'j': 0, 'R': (1,), 't': 1.}])
        chain.set_onsite({'a': (0.5, -0.5)})
        en = np.linalg.eigvalsh(chain.get_ham((0.3,)))
        self.assertFalse(np.isclose(en[0], en[1]))
        self.assertTrue(np.isclose(en[1] - en[0], 1.))

    def test_hermitian(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.), (0., 1.)])
        sq = KSpace(lat, spin=True)
        sq.set_hopping([{'i': 0, 'j': 0, 'R': (1, 0), 't': PAULI['0'] + 1j*0.3*PAULI['y']},
                              {'i': 0, 'j': 0, 'R': (0, 1), 't': PAULI['0'] - 1j*0.3*PAULI['x']}])
        for k in ([0.3, 0.4], [0., 0.], [np.pi, 0.]):
            ham = sq.get_ham(k)
            self.assertTrue(np.allclose(ham, ham.conj().T))

    def test_rashba_splits_away_from_trim_points(self):
        lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.), (0., 1.)])
        sq = KSpace(lat, spin=True)
        sq.set_hopping([{'i': 0, 'j': 0, 'R': (1, 0), 't': PAULI['0'] + 1j*0.3*PAULI['y']},
                              {'i': 0, 'j': 0, 'R': (0, 1), 't': PAULI['0'] - 1j*0.3*PAULI['x']}])
        # time-reversal-invariant momenta: Kramers degeneracy must survive.
        for k in ([0., 0.], [np.pi, 0.], [0., np.pi], [np.pi, np.pi]):
            en = np.linalg.eigvalsh(sq.get_ham(k))
            self.assertTrue(np.isclose(en[0], en[1]))
        # generic k-point: Rashba lifts the spin degeneracy.
        en = np.linalg.eigvalsh(sq.get_ham([0.3, 0.4]))
        self.assertFalse(np.isclose(en[0], en[1]))

    def test_kane_mele_matches_decoupled_spin_sectors(self):
        # Intrinsic (sz-conserving) spin-orbit coupling on honeycomb: the
        # spinful model must be exactly equivalent to two decoupled
        # spinless models with opposite next-nearest-neighbor chirality.
        t1, lam = 1., 0.1
        km, up, down = kane_mele(t1, lam)
        for k in ([0.3, -0.7], [1.1, 0.4], [0., 0.]):
            en_full = np.sort(np.linalg.eigvalsh(km.get_ham(k)))
            en_sectors = np.sort(np.concatenate([np.linalg.eigvalsh(up.get_ham(k)),
                                                                       np.linalg.eigvalsh(down.get_ham(k))]))
            self.assertTrue(np.allclose(en_full, en_sectors))
            # time-reversal symmetry -> Kramers degeneracy at every k.
            self.assertTrue(np.isclose(en_full[0], en_full[1]))
            self.assertTrue(np.isclose(en_full[2], en_full[3]))


def haldane_sector(t1, lam, sign):
    lat = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}],
                          prim_vec=[(2*DX, 0.), (DX, 1.5)])
    s = KSpace(lat)
    s.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t1},
                          {'i': 0, 'j': 1, 'R': (-1, 0), 't': t1},
                          {'i': 0, 'j': 1, 'R': (0, -1), 't': t1}])
    for R in [(1, 0), (0, 1), (1, -1)]:
        s.set_hopping([{'i': 0, 'j': 0, 'R': R, 't': 1j*sign*lam}])
        s.set_hopping([{'i': 1, 'j': 1, 'R': R, 't': -1j*sign*lam}])
    return s


def kane_mele(t1, lam):
    unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
    prim_vec = [(2*DX, 0.), (DX, 1.5)]
    lat = lattice(unit_cell=unit_cell, prim_vec=prim_vec)
    km = KSpace(lat, spin=True)
    km.set_hopping([{'i': 0, 'j': 1, 'R': (0, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (-1, 0), 't': t1},
                            {'i': 0, 'j': 1, 'R': (0, -1), 't': t1}])
    for R in [(1, 0), (0, 1), (1, -1)]:
        km.set_hopping([{'i': 0, 'j': 0, 'R': R, 't': 1j*lam*PAULI['z']}])
        km.set_hopping([{'i': 1, 'j': 1, 'R': R, 't': -1j*lam*PAULI['z']}])
    return km, haldane_sector(t1, lam, +1), haldane_sector(t1, lam, -1)


GRAPHENE_LIST_HOP = [{'i': 0, 'j': 1, 'R': (0, 0), 't': 1.},
                                    {'i': 0, 'j': 1, 'R': (-1, 0), 't': 1.},
                                    {'i': 0, 'j': 1, 'R': (0, -1), 't': 1.}]


class TestRibbon(unittest.TestCase):
    '''
    Unittest of the free function **ribbon**.
    '''
    def graphene_lat(self):
        unit_cell = [{'tag': 'a', 'r0': (0., 0.)}, {'tag': 'b', 'r0': (DX, DY)}]
        prim_vec = [(2*DX, 0.), (DX, 1.5)]
        return lattice(unit_cell=unit_cell, prim_vec=prim_vec)

    def test_errors(self):
        lat = self.graphene_lat()
        self.assertRaises(TypeError, ribbon, 0, GRAPHENE_LIST_HOP, 10)
        self.assertRaises(ValueError, ribbon, lat, GRAPHENE_LIST_HOP, -1)
        self.assertRaises(ValueError, ribbon, lat, GRAPHENE_LIST_HOP, 10, direction=2)
        self.assertRaises(TypeError, ribbon, lat, GRAPHENE_LIST_HOP, 10, direction=1.5)
        lat_1d = lattice(unit_cell=[{'tag': 'a', 'r0': (0., 0.)}], prim_vec=[(1., 0.)])
        self.assertRaises(ValueError, ribbon, lat_1d, [], 10)

    def test_shape(self):
        lat = self.graphene_lat()
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=10)
        self.assertEqual(rib.dim, 1)
        self.assertEqual(rib.n_sites, 20)
        self.assertEqual(rib.norb, 20)

    def test_hermitian_and_bulk_bandwidth(self):
        lat = self.graphene_lat()
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=10)
        for k in (0., 1.3, np.pi):
            ham = rib.get_ham((k,))
            self.assertTrue(np.allclose(ham, ham.conj().T))
        en = rib.mesh_bands(nk=200)
        # a wide-enough ribbon reproduces the bulk graphene bandwidth 3t.
        self.assertTrue(np.isclose(en.max(), 3., atol=0.05))
        self.assertTrue(np.isclose(en.min(), -3., atol=0.05))

    def test_zigzag_edge_states(self):
        # the hallmark zigzag-graphene-ribbon result: a partially-flat band
        # pinned at E=0 over part of the 1D Brillouin zone.
        lat = self.graphene_lat()
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=30, direction=1)
        en = rib.mesh_bands(nk=300)
        n_zero_modes = np.sum(np.any(np.abs(en) < 1e-6, axis=1))
        self.assertTrue(n_zero_modes > 50)

    def test_onsite_applied_per_row(self):
        lat = self.graphene_lat()
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=3, onsite={'a': 1., 'b': -1.})
        self.assertTrue(np.allclose(rib.onsite[0::2], 1.))
        self.assertTrue(np.allclose(rib.onsite[1::2], -1.))

    def test_spin(self):
        lat = self.graphene_lat()
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=5, spin=True)
        self.assertTrue(rib.spin)
        self.assertEqual(rib.norb, 20)
        en = np.linalg.eigvalsh(rib.get_ham((0.3,)))
        self.assertEqual(len(en) % 2, 0)
        # spin-independent hopping -> every level doubly degenerate.
        self.assertTrue(np.allclose(en[0::2], en[1::2]))

    def test_matches_large_real_space_flake(self):
        lat = self.graphene_lat()
        width = 8
        rib = ribbon(lat, GRAPHENE_LIST_HOP, width=width, direction=1)
        n1 = 40
        lat2 = self.graphene_lat()
        lat2.get_lattice(n1=n1, n2=width)
        sys = system(lat=lat2)
        sys.set_hopping([{'n': 1, 't': 1.}])
        sys.get_ham()
        sys.get_eig()
        flake_en = sys.en.real
        ribbon_en = rib.mesh_bands(nk=n1).ravel()
        for tol in (0.05, 0.1):
            n_flake = np.sum(np.abs(flake_en) < tol)
            n_ribbon = np.sum(np.abs(ribbon_en) < tol)
            self.assertTrue(abs(n_flake - n_ribbon) <= 4)


if __name__ == '__main__':
    unittest.main()
