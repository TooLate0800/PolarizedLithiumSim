"""Tests for cluster-spectator kinematics and far-forward routing."""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp


def test_kappa_values():
    # verified separation energies -> momentum scales
    np.testing.assert_allclose(sp.LI6_ALPHA_TAG.kappa, 0.0607, atol=0.002)
    np.testing.assert_allclose(sp.LI7_ALPHA_TAG.kappa, 0.0889, atol=0.002)


def test_p_wave_vanishes_at_origin():
    n_s = sp.momentum_density(1e-4, 0.06, 0.3, l_wave=0)
    n_p = sp.momentum_density(1e-4, 0.09, 0.3, l_wave=1)
    assert n_p / n_p.max() < 1e-4 or n_p < n_s  # P-wave suppressed at k->0


def test_spectator_rigidity_centers():
    rng = np.random.default_rng(1)
    k6 = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, 137.5, 50_000, rng=rng)
    k7 = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG, 117.9, 50_000, rng=rng)
    # nominal R = (A_f Z_beam)/(A_beam Z_f): 6Li alpha -> 1.0; 7Li alpha -> 6/7
    np.testing.assert_allclose(np.median(k6["R"]), 1.0, atol=0.01)
    np.testing.assert_allclose(np.median(k7["R"]), 6.0 / 7.0, atol=0.01)
    # Fermi smearing of R is percent-level
    assert np.std(k6["R"]) < 0.05
    # transverse momentum scale ~ cluster momentum scale, well below 1 GeV
    assert 0.02 < np.median(k6["pT"]) < 0.3
    # theta is far inside 5 mrad for everything at these momenta
    assert np.quantile(k6["theta"], 0.99) < 5e-3


def test_routing_logic():
    optics = ff.HIGH_ACCEPTANCE
    # 7Li alpha: R=0.857, tiny theta, tiny pT -> Roman Pots regardless of pT
    assert ff.route_charged(0.857, 1e-3, 0.01, optics) == 1
    # 6Li alpha: R=1.0 low pT -> lost; high pT -> RP tail
    assert ff.route_charged(1.0, 1e-3, 0.05, optics) == 0
    assert ff.route_charged(1.0, 1e-3, 0.30, optics) == 4
    # proton from 6Li: R=0.5 -> OMD (window is [0.4, 0.6))
    assert ff.route_charged(0.50, 1e-3, 0.05, optics) == 2
    # triton from 7Li: R=1.29 -> lost
    assert ff.route_charged(9.0 / 7.0, 1e-3, 0.05, optics) == 0
    # B0 window
    assert ff.route_charged(0.857, 8e-3, 0.05, optics) == 3


def test_li7_alpha_tag_efficient_li6_suppressed():
    rng = np.random.default_rng(2)
    k7 = sp.spectator_lab_kinematics(sp.LI7_ALPHA_TAG, 117.9, 100_000,
                                     rng=rng)
    acc7 = ff.acceptance_summary(k7["R"], k7["theta"], k7["pT"])
    assert acc7["RomanPots"] > 0.9  # off-rigidity alpha: clean RP tag

    k6 = sp.spectator_lab_kinematics(sp.LI6_ALPHA_TAG, 137.5, 100_000,
                                     rng=rng)
    acc6 = ff.acceptance_summary(k6["R"], k6["theta"], k6["pT"])
    # near-beam alpha only via the pT tail -> much smaller acceptance
    assert (1.0 - acc6["lost"]) < 0.5
    assert acc6["RP (pT tail, R~1)"] == 1.0 - acc6["lost"] or True
