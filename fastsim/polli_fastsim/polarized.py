"""Polarized / tensor / gluonometric inputs (TOY models + hooks).

Every quantity here is a *scenario input*, not a prediction. The toy shapes
have roughly the right magnitude and x-dependence so that FOM machinery can
be exercised end-to-end; Phase-1 steps 1.3-1.5 replace them with:
  - g1p, g1n: JAM / DSSV grids (LHAPDF)
  - polarized EMC ratio: digitized Cloet-Bentz-Thomas curves (PLB 642:210)
  - b1: Cosyn-Dong-Kumano-Song deuteron convolution + 6Li cluster model
  - Delta(x,Q2): lattice-moment-normalized model (Detmold-Shanahan) or
    scenario values Delta/F1 in {1e-3 ... 1e-2}
"""

import numpy as np

from .structure import ToyF2, r_sigma_lt


class ToyG1:
    """g1 = A1 * F1 with toy A1(x) shapes (TOY; replace with JAM/DSSV)."""

    def __init__(self, base=None):
        self.base = base or ToyF2()

    def a1p(self, x):
        # x^0.7 tracks moderate/high-x world data; saturates below 1
        return np.clip(np.power(np.asarray(x, dtype=float), 0.7), 0.0, 1.0)

    def a1n(self, x):
        x = np.asarray(x, dtype=float)
        # small negative at low/mid x, positive rise at high x
        return -0.07 * np.power(1.0 - x, 2.0) + 0.8 * np.power(x, 2.2)

    def _f1(self, f2, x, q2):
        return f2 / (2.0 * x * (1.0 + r_sigma_lt(x, q2)))

    def g1p(self, x, q2):
        return self.a1p(x) * self._f1(self.base.f2p(x, q2), x, q2)

    def g1n(self, x, q2):
        return self.a1n(x) * self._f1(self.base.f2n(x, q2), x, q2)

    def g1_nucleus(self, ion, x, q2, medium_ratio=None):
        """g1A = P_p g1p + P_n g1n (naive), times optional medium ratio.

        `medium_ratio(x)` is the polarized-EMC modification DR(x); the
        observable of interest is precisely deviations of DR from 1.
        """
        g1 = ion.eff_pol_p * self.g1p(x, q2) + ion.eff_pol_n * self.g1n(x, q2)
        if medium_ratio is not None:
            g1 = g1 * medium_ratio(x)
        return g1


def toy_polarized_emc_ratio(x):
    """TOY stand-in for the Cloet-Bentz-Thomas polarized-EMC ratio DR(x).

    CBT predict the polarized EMC effect to be ~2x the unpolarized one in
    the valence region (ratio dipping to ~0.8 around x~0.6-0.7 for 7Li).
    Shape below is qualitative only -- digitize the published curve for
    real projections.
    """
    x = np.asarray(x, dtype=float)
    return 1.0 - 0.25 * np.power(x, 1.5) * np.exp(-2.0 * (x - 0.75) ** 2)


def toy_b1(x, q2, f1):
    """TOY b1 with HERMES-like magnitude: ~ +1e-2*F1 at x~0.01 crossing
    zero near x~0.3-0.45 (cf. HERMES PRL 95:242001). Returns b1 (absolute).
    """
    x = np.asarray(x, dtype=float)
    shape = 0.01 * np.power(np.maximum(x, 1e-6), -0.2) * (1.0 - x / 0.35)
    return shape * np.exp(-3.0 * x) * f1


def toy_delta_gluon(x, q2, f1, scale=1e-3):
    """TOY double-helicity-flip Delta(x,Q2) = scale * F1 * (1-x)^4 * x^0.3.

    Delta is alpha_s-suppressed and unknown; `scale` sets Delta/F1 at its
    maximum. Use scale in {1e-3, 3e-3, 1e-2} as discovery scenarios, to be
    replaced by lattice-moment-anchored models.
    """
    x = np.asarray(x, dtype=float)
    return scale * f1 * np.power(x, 0.3) * np.power(1.0 - x, 4.0)
