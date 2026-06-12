"""Toy-MC closure of the analytic statistical-error formulas (step 1.3.4).

Validates that the per-bin error estimators in asymmetries.py reproduce
the spread of explicit Poisson/sampling pseudo-experiments.
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim.asymmetries import err_azz, err_cos2phi_amplitude


def test_azz_estimator_closure():
    rng = np.random.default_rng(42)
    n_tot, eps, pzz = 1_000_000, 0.002, 0.8  # physics Azz = 2*eps
    trials = 400
    c = np.array([1.0, -2.0, 1.0])
    lam = n_tot / 3.0 * (1.0 + c[None, :] * eps * pzz)
    counts = rng.poisson(np.broadcast_to(lam, (trials, 3)))
    est = ((counts[:, 0] + counts[:, 2] - 2 * counts[:, 1])
           / counts.sum(axis=1)) / pzz
    # mean recovers the physics asymmetry 2*eps
    np.testing.assert_allclose(est.mean(), 2 * eps,
                               atol=4 * est.std() / np.sqrt(trials))
    # spread matches the analytic formula sqrt(2/N)/Pzz within 10%
    expected = err_azz(n_tot, pzz)
    assert 0.9 < est.std() / expected < 1.1


def test_cos2phi_moment_estimator_closure():
    rng = np.random.default_rng(7)
    n, amp, pzz = 100_000, 0.01, 1.0
    trials = 300
    stds, means = [], []
    for _ in range(trials):
        # accept-reject sampling of (1 + amp*cos 2phi)/2pi
        phi = rng.uniform(0, 2 * np.pi, size=int(n * 1.3))
        keep = rng.uniform(0, 1 + amp, size=phi.size) < 1 + amp * np.cos(2 * phi)
        phi = phi[keep][:n]
        means.append(2.0 * np.mean(np.cos(2.0 * phi)))
    means = np.array(means)
    np.testing.assert_allclose(means.mean(), amp,
                               atol=4 * means.std() / np.sqrt(trials))
    expected = err_cos2phi_amplitude(n, pzz)
    assert 0.9 < means.std() / expected < 1.1
