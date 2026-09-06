"""Independent thermodynamic and logical-channel checks for the candidate."""
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm

from rh1_passive_reference import Reference, gaussian_nominal_error


class PassiveReferenceTests(unittest.TestCase):
    def setUp(self):
        self.ref = Reference(k0=1.1, temperature=.3)

    def test_force_is_actual_gibbs_expectation(self):
        x = np.array([[0., 1.], [1., 0.]])
        z = np.diag([1., -1.])
        for s, delta in ((-.3, .2), (.1, -.4), (.2, .1)):
            h = (x + (delta-s)*z)/2
            rho = expm(-h/self.ref.temperature)
            rho /= np.trace(rho)
            force = -self.ref.k0*s + np.trace(rho@z).real/2
            self.assertAlmostEqual(self.ref.force(s, delta), force, places=13)
            step = 1e-5
            numerical = -(self.ref.potential(s+step, delta)
                          - self.ref.potential(s-step, delta))/(2*step)
            self.assertAlmostEqual(numerical, force, places=9)

    def test_stiffness_and_rest_from_both_sides(self):
        for delta in (-.2, 0., .2):
            root = self.ref.rest(delta)
            self.assertLess(root*delta, 1e-15)
            for displacement in (-.03, .03):
                self.assertLess(self.ref.force(root+displacement, delta)*displacement, 0)
            step = 1e-5
            slope = -(self.ref.force(root+step, delta)
                      - self.ref.force(root-step, delta))/(2*step)
            self.assertAlmostEqual(slope, self.ref.stiffness(root, delta), places=9)

    def test_partition_function_produces_local_variance(self):
        # A direct non-Gaussian canonical integral checks the harmonic limit.
        n, temperature = 10000, .3
        ref = Reference(k0=1.1, temperature=temperature, number=n)
        sigma = np.sqrt(temperature/ref.stiffness(0., 0.))
        baseline = ref.potential(0., 0.)
        density = lambda y: np.exp(-(ref.potential(sigma*y, 0.)-baseline)/temperature)
        norm = quad(density, -10, 10, epsabs=1e-11)[0]
        variance = sigma**2*quad(lambda y:y*y*density(y), -10, 10, epsabs=1e-11)[0]/norm
        self.assertLess(abs(variance/sigma**2-1), 2e-4)

    def test_gaussian_error_includes_coherence_loss(self):
        mean, variance = .012, .0004
        actual = quad(lambda y: (2/3)*np.sin(2*(mean+np.sqrt(variance)*y))**2
                      * np.exp(-y*y/2)/np.sqrt(2*np.pi), -12, 12)[0]
        self.assertAlmostEqual(actual, gaussian_nominal_error(mean, variance), places=13)
        self.assertGreater(gaussian_nominal_error(0., variance), 0.)


if __name__ == '__main__':
    unittest.main()
