"""DIS kinematics for e+ion collisions (per-nucleon variables).

Conventions (EIC): ion beam along +z, electron beam along -z, massless
approximation throughout. theta is the polar angle from +z (ion direction),
so the scattered electron sits at theta ~ pi (negative eta).

All structure-function kinematics use the electron-NUCLEON system:
s = 4 E_e E_N with E_N the ion energy per nucleon; x = x_Bjorken per
nucleon in [0, ~A] but we restrict to x <= 1 for this study.
"""

import numpy as np


def sqrt_s(electron_energy, ion_energy_per_nucleon):
    return np.sqrt(4.0 * electron_energy * ion_energy_per_nucleon)


def q2_from_xy(x, y, s):
    return s * x * y


def y_from_xq2(x, q2, s):
    return q2 / (s * x)


def w2(x, q2, m_n=0.9383):
    """Invariant mass^2 of the hadronic system (per-nucleon target)."""
    return m_n**2 + q2 * (1.0 - x) / np.maximum(x, 1e-12)


def scattered_electron(x, y, s, electron_energy):
    """Lab-frame scattered electron (E', theta, eta).

    From k=(Ee,0,0,-Ee), P=(EN,0,0,EN):
      Q2 = 2 Ee E' (1 + cos(theta)),   y = 1 - (E'/2Ee)(1 - cos(theta))
    =>  E'(1+cos) = Q2/(2Ee),          E'(1-cos) = 2 Ee (1-y)
    """
    q2 = q2_from_xy(x, y, s)
    a = q2 / (2.0 * electron_energy)       # E'(1+cos)
    b = 2.0 * electron_energy * (1.0 - y)  # E'(1-cos)
    e_prime = 0.5 * (a + b)
    cos_theta = np.where(e_prime > 0, (a - b) / (2.0 * e_prime), -1.0)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.arccos(cos_theta)
    # eta = -ln tan(theta/2); guard the exactly-backward direction
    tan_half = np.tan(np.minimum(theta, np.pi - 1e-9) / 2.0)
    eta = -np.log(np.maximum(tan_half, 1e-12))
    return e_prime, theta, eta


def kinematic_mask(
    x,
    q2,
    s,
    q2_min=1.0,
    y_min=0.01,
    y_max=0.95,
    w2_min=10.0,
    x_max=1.0,
):
    """Standard inclusive-DIS phase-space cuts (Yellow-Report-like)."""
    y = y_from_xq2(x, q2, s)
    return (
        (q2 >= q2_min)
        & (y >= y_min)
        & (y <= y_max)
        & (w2(x, q2) >= w2_min)
        & (x <= x_max)
    )


def log_grid(x_range=(1e-4, 1.0), q2_range=(1.0, 1e3), nx=40, nq2=30):
    """Bin edges and centers, log-spaced, for (x, Q2) maps."""
    x_edges = np.logspace(np.log10(x_range[0]), np.log10(x_range[1]), nx + 1)
    q2_edges = np.logspace(np.log10(q2_range[0]), np.log10(q2_range[1]), nq2 + 1)
    x_c = np.sqrt(x_edges[:-1] * x_edges[1:])
    q2_c = np.sqrt(q2_edges[:-1] * q2_edges[1:])
    return x_edges, q2_edges, x_c, q2_c
