"""Unpolarized structure functions: toy parameterizations + nuclear builder.

The TOY F2 below is a two-component (valence + sea) form, anchored by eye
to world ep data: F2p(1e-3, 10) ~ 1.2, F2p(0.1, 10) ~ 0.36, F2p(0.4, 10) ~ 0.12.
It is adequate for phase-space maps and ~factor-1.5 rate estimates ONLY.

Step 1.3 of the Phase-1 plan replaces this with LHAPDF (CT18 + EPPS21 or
nNNPDF) evaluated inside eic-shell, or the `parton` pip package locally.
Keep the F2Source interface so the swap is one line in the scripts.
"""

import numpy as np

ALPHA_EM = 1.0 / 137.036
GEV2_TO_PB = 0.3894e9  # (hbar c)^2 in GeV^2 * pb


class ToyF2:
    """Crude but smooth F2p / F2n with mild log Q2 evolution (TOY)."""

    def f2p(self, x, q2):
        x = np.asarray(x, dtype=float)
        q2 = np.asarray(q2, dtype=float)
        lq = np.log(np.maximum(q2, 1.1) / 0.04)  # ln(Q2/Lambda^2), Lambda=0.2
        # sea: low-x rise steepens slowly with Q2 (HERA-like lambda ~ 0.05 lnQ2)
        lam = 0.045 * lq
        sea = 0.20 * np.power(x, -lam) * np.power(1.0 - x, 7.0)
        val = 1.05 * np.power(x, 0.55) * np.power(1.0 - x, 3.0)
        return sea + val

    def f2n_over_f2p(self, x):
        # simple fit to NMC-like d/p ratio behavior
        return np.clip(1.0 - 0.75 * np.asarray(x, dtype=float), 0.25, 1.0)

    def f2n(self, x, q2):
        return self.f2p(x, q2) * self.f2n_over_f2p(x)


def r_sigma_lt(x, q2):
    """R = sigma_L/sigma_T, simplified R1990-like magnitude."""
    q2 = np.asarray(q2, dtype=float)
    return 0.18 / (1.0 + q2 / 50.0)


class NuclearF2:
    """Whole-nucleus F2A = Z*F2p + N*F2n, with an optional EMC-ratio hook.

    Per-nucleon F2 is F2A/A. The emc_ratio callable r(x) multiplies the
    isoscalar combination; default None means r = 1 (no medium effect),
    fine for rates (the EMC effect is a <20% shape effect here).
    """

    def __init__(self, ion, base=None, emc_ratio=None):
        self.ion = ion
        self.base = base or ToyF2()
        self.emc_ratio = emc_ratio

    def f2a(self, x, q2):
        f2 = self.ion.Z * self.base.f2p(x, q2) + self.ion.N * self.base.f2n(x, q2)
        if self.emc_ratio is not None:
            f2 = f2 * self.emc_ratio(x)
        return f2

    def f1a(self, x, q2):
        """F1 via the Callan-Gross relation modified by R (massless)."""
        r = r_sigma_lt(x, q2)
        return self.f2a(x, q2) / (2.0 * np.asarray(x, dtype=float) * (1.0 + r))


def dsigma_dx_dq2(x, q2, s, f2, fl=None):
    """NC DIS double-differential cross section [pb/GeV^2].

    d2sigma/dxdQ2 = 4 pi alpha^2 / (x Q^4) [ (1 - y + y^2/2) F2 - y^2/2 FL ]
    FL defaults to F2 * R/(1+R).
    """
    x = np.asarray(x, dtype=float)
    q2 = np.asarray(q2, dtype=float)
    y = q2 / (s * x)
    if fl is None:
        r = r_sigma_lt(x, q2)
        fl = f2 * r / (1.0 + r)
    bracket = (1.0 - y + 0.5 * y * y) * f2 - 0.5 * y * y * fl
    xsec = 4.0 * np.pi * ALPHA_EM**2 / (x * q2**2) * bracket
    return np.maximum(xsec, 0.0) * GEV2_TO_PB
