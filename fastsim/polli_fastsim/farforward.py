"""Far-forward detector windows and fragment routing (Phase-1 stand-in).

Windows fetch-verified 2026-06-12 (YR detector matrix; arXiv:2108.08314
Table I; arXiv:2409.02811; details in plans/03 step 2.2):

  Roman Pots   z = 26/28 m   theta < 5 mrad   R in [0.60, 0.95]
               near-beam (|R-1| < 0.05): only pT > pT_cut (10-sigma optics)
  OMD          z = 22.5/24.5 m  theta < 5 mrad  R in [0.40, 0.60], no pT cut
  B0           5.5 < theta < 20 mrad, any charged
  ZDC          theta < 4 mrad, neutrals
  no coverage  R > ~1.05 (bends less than beam), or 5-5.5 mrad gap

pT_cut: 0.20 GeV (high-acceptance optics) / 0.45 GeV (high-divergence).
ASSUMPTION (flagged in plans/04 #11): off-rigidity tracks inside the RP
window are dispersion-separated from the beam, so no pT cut applies to
them; only Phase-2 optics can refine this.
"""

from dataclasses import dataclass

import numpy as np

THETA_RP_MAX = 5.0e-3
THETA_B0_MIN, THETA_B0_MAX = 5.5e-3, 20.0e-3
THETA_ZDC_MAX = 4.0e-3
RP_R_WINDOW = (0.60, 0.95)
OMD_R_WINDOW = (0.40, 0.60)
NEAR_BEAM_BAND = 0.05  # |R-1| below this: inside the beam envelope unless pT


@dataclass(frozen=True)
class Optics:
    name: str
    pt_cut_near_beam: float  # GeV


HIGH_ACCEPTANCE = Optics("high-acceptance", 0.20)
HIGH_DIVERGENCE = Optics("high-divergence", 0.45)


def route_charged(R, theta, pT, optics=HIGH_ACCEPTANCE):
    """Classify charged fragments into far-forward systems.

    Returns an integer array: 0 lost, 1 RP, 2 OMD, 3 B0, 4 RP-near-beam
    (R ~ 1, accepted only above the optics pT cut).
    """
    R = np.asarray(R, dtype=float)
    theta = np.asarray(theta, dtype=float)
    pT = np.asarray(pT, dtype=float)
    out = np.zeros(R.shape, dtype=int)

    in_b0 = (theta >= THETA_B0_MIN) & (theta <= THETA_B0_MAX)
    out[in_b0] = 3

    small = theta < THETA_RP_MAX
    near = np.abs(R - 1.0) < NEAR_BEAM_BAND
    in_rp = small & ~near & (R >= RP_R_WINDOW[0]) & (R <= RP_R_WINDOW[1])
    in_omd = small & (R >= OMD_R_WINDOW[0]) & (R < OMD_R_WINDOW[1])
    rp_tail = small & near & (pT > optics.pt_cut_near_beam)
    out[in_omd] = 2
    out[in_rp] = 1
    out[rp_tail] = 4
    return out


ROUTE_LABELS = {0: "lost", 1: "RomanPots", 2: "OMD", 3: "B0",
                4: "RP (pT tail, R~1)"}


def acceptance_summary(R, theta, pT, optics=HIGH_ACCEPTANCE):
    """Fraction of spectators in each far-forward system."""
    route = route_charged(R, theta, pT, optics)
    n = float(len(route))
    return {label: float(np.sum(route == code)) / n
            for code, label in ROUTE_LABELS.items()}
