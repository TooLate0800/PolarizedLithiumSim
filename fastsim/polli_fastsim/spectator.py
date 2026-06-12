"""Cluster-spectator kinematics for e+6Li and e+7Li tagging.

Model: the lithium ground state is a two-cluster system (6Li = alpha + d,
S-wave dominant; 7Li = alpha + t, P-wave, 3/2-). DIS strikes one cluster;
the other is a spectator carrying the internal relative momentum k. We
sample k from simple analytic momentum densities, boost to the lab, and
report the far-forward routing variables (pT, theta, rigidity ratio R).

Momentum-density models (two-parameter, crude BY DESIGN — the high-k tail
is the dominant model uncertainty and exactly what VMC densities/BeAGLE
must eventually pin down; see plans/02 step 1.5.3):

  S-wave (6Li):  psi(k) ∝ 1/(k^2+kappa^2) - 1/(k^2+beta^2)   (Hulthen form)
  P-wave (7Li):  psi(k) ∝ k / [(k^2+kappa^2)(k^2+beta^2)]    (vanishes at 0)

kappa = sqrt(2*mu*S_alpha) from the verified separation energies:
  6Li -> alpha+d: S = 1.474 MeV  => kappa = 60.7 MeV
  7Li -> alpha+t: S = 2.467 MeV  => kappa = 88.9 MeV
beta is the short-range scale (default 0.30 GeV, scan 0.20-0.40 to span
the tail uncertainty; deuteron analogy: beta/kappa ~ 5.7).

Upgrade path: replace n(k) with VMC two-cluster overlap densities
(R.B. Wiringa et al., tables available from ANL) — same interface.
"""

from dataclasses import dataclass

import numpy as np

M_U = 0.93149  # amu [GeV]
MASSES = {"p": 0.93827, "n": 0.93957, "d": 1.87561, "t": 2.80892,
          "3He": 2.80839, "alpha": 3.72738}


@dataclass(frozen=True)
class ClusterChannel:
    """Beam nucleus -> struck cluster + spectator cluster."""
    name: str            # e.g. "6Li -> d* + alpha(spec)"
    beam_A: int
    beam_Z: int
    spectator: str       # key into MASSES
    spectator_A: int
    spectator_Z: int
    separation_energy: float  # [GeV]
    l_wave: int          # 0 (S) or 1 (P) relative motion

    @property
    def m_spec(self):
        return MASSES[self.spectator]

    @property
    def kappa(self):
        """Bound-state momentum scale sqrt(2 mu S) [GeV]."""
        m_beam = self.beam_A * M_U  # adequate for mu at this precision
        m_other = m_beam - self.m_spec
        mu = self.m_spec * m_other / (self.m_spec + m_other)
        return np.sqrt(2.0 * mu * self.separation_energy)


LI6_ALPHA_TAG = ClusterChannel(
    "6Li: DIS on d-cluster, alpha spectator", 6, 3, "alpha", 4, 2,
    1.4743e-3, l_wave=0)
LI6_D_TAG = ClusterChannel(
    "6Li: DIS on alpha-cluster, d spectator", 6, 3, "d", 2, 1,
    1.4743e-3, l_wave=0)
LI7_ALPHA_TAG = ClusterChannel(
    "7Li: DIS on t-cluster, alpha spectator", 7, 3, "alpha", 4, 2,
    2.4670e-3, l_wave=1)
LI7_T_TAG = ClusterChannel(
    "7Li: DIS on alpha-cluster, triton spectator", 7, 3, "t", 3, 1,
    2.4670e-3, l_wave=1)

CHANNELS = (LI6_ALPHA_TAG, LI6_D_TAG, LI7_ALPHA_TAG, LI7_T_TAG)


def momentum_density(k, kappa, beta, l_wave):
    """Unnormalized n(k) (NOT including the k^2 phase-space factor)."""
    k = np.asarray(k, dtype=float)
    if l_wave == 0:
        psi = 1.0 / (k * k + kappa * kappa) - 1.0 / (k * k + beta * beta)
    else:
        psi = k / ((k * k + kappa * kappa) * (k * k + beta * beta))
    return psi * psi


def sample_k(channel, n, beta=0.30, k_max=1.5, rng=None):
    """Sample |k| from k^2 * n(k) on [0, k_max] GeV by inverse-CDF on a grid,
    plus isotropic angles. Returns (kx, ky, kz) arrays."""
    rng = rng or np.random.default_rng(20260612)
    grid = np.linspace(1e-4, k_max, 30000)
    pdf = grid * grid * momentum_density(grid, channel.kappa, beta,
                                         channel.l_wave)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    k = np.interp(rng.uniform(size=n), cdf, grid)
    cos_t = rng.uniform(-1.0, 1.0, size=n)
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n)
    sin_t = np.sqrt(1.0 - cos_t**2)
    return k * sin_t * np.cos(phi), k * sin_t * np.sin(phi), k * cos_t


def spectator_lab_kinematics(channel, p_per_nucleon, n=200_000, beta=0.30,
                             rng=None):
    """Boost spectator (rest-frame momentum -k of the struck cluster) to the
    lab. Returns dict of arrays: pT, theta [rad], p_lab, R (rigidity ratio
    vs beam), xL (= p_lab / (A_spec * p_per_nucleon))."""
    kx, ky, kz = sample_k(channel, n, beta=beta, rng=rng)
    m = channel.m_spec
    e_rest = np.sqrt(m * m + kx * kx + ky * ky + kz * kz)
    # beam boost (per-nucleon momentum sets the velocity of the nucleus)
    m_beam = channel.beam_A * M_U
    p_beam = channel.beam_A * p_per_nucleon
    e_beam = np.sqrt(p_beam**2 + m_beam**2)
    gamma = e_beam / m_beam
    gbeta = p_beam / m_beam
    pz_lab = gamma * kz + gbeta * e_rest
    pt = np.sqrt(kx * kx + ky * ky)
    p_lab = np.sqrt(pt * pt + pz_lab * pz_lab)
    theta = np.arctan2(pt, pz_lab)
    rigidity_beam = p_beam / channel.beam_Z
    rigidity_spec = p_lab / channel.spectator_Z
    return {
        "pT": pt,
        "theta": theta,
        "p_lab": p_lab,
        "R": rigidity_spec / rigidity_beam,
        "xL": p_lab / (channel.spectator_A * p_per_nucleon),
        "k": np.sqrt(kx**2 + ky**2 + kz**2),
    }
