import unittest
import numpy as np

import tbee.dos as dos


class TestDos(unittest.TestCase):
    '''
    Unittest of module **dos**.
    '''
    def test_errors(self):
        self.assertRaises(ValueError, dos.density_of_states, np.array([]))
        self.assertRaises(TypeError, dos.density_of_states, np.array([1.]), broadening='a')
        self.assertRaises(ValueError, dos.density_of_states, np.array([1.]), broadening=-1.)
        self.assertRaises(TypeError, dos.density_of_states, np.array([1.]), kernel=0)
        self.assertRaises(ValueError, dos.density_of_states, np.array([1.]), kernel='z')
        self.assertRaises(ValueError, dos.density_of_states, np.array([1.]), e_grid=np.array([]))

    def test_normalization_gaussian(self):
        rng = np.random.default_rng(0)
        energies = rng.uniform(-2., 2., 500)
        e_grid, rho = dos.density_of_states(energies, broadening=0.1)
        self.assertTrue(np.isclose(np.trapezoid(rho, e_grid), len(energies), rtol=1e-3))

    def test_normalization_lorentzian(self):
        rng = np.random.default_rng(1)
        energies = rng.uniform(-2., 2., 500)
        e_grid, rho = dos.density_of_states(energies, broadening=0.1, kernel='lorentzian',
                                                                  e_grid=np.linspace(-20., 20., 4001))
        self.assertTrue(np.isclose(np.trapezoid(rho, e_grid), len(energies), rtol=1e-2))

    def test_single_level_peak(self):
        e_grid, rho = dos.density_of_states(np.array([0.]), broadening=0.2,
                                                                  e_grid=np.linspace(-2., 2., 401))
        self.assertEqual(e_grid[np.argmax(rho)], 0.)

    def test_complex_input_uses_real_part(self):
        e_grid1, rho1 = dos.density_of_states(np.array([1., -1.]), broadening=0.2)
        e_grid2, rho2 = dos.density_of_states(np.array([1.+5j, -1.-5j]), broadening=0.2)
        self.assertTrue(np.allclose(e_grid1, e_grid2))
        self.assertTrue(np.allclose(rho1, rho2))


if __name__ == '__main__':
    unittest.main()
