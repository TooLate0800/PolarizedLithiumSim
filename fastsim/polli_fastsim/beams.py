"""Beam species and EIC energy configurations.

Ion top momenta scale with magnetic rigidity from the 275 GeV proton ring:
    p/nucleon ~= 275 GeV * (Z/A)
This reproduces d (137.5, CDR rounds to "135") and Au (110 GeV/u). 3He is
183 GeV/u in the current design (CDR Sec. 5.5, 6-snake spin preservation);
the 166 GeV/u seen in older papers is the eRHIC-era 250-GeV-proton legacy.
Li values below match the EPIOS white paper (Atoian et al.,
arXiv:2510.10794, PRC 113:060501: ~138 GeV/u 6Li, ~117 GeV/u 7Li);
final numbers still need C-AD blessing.
"""

from dataclasses import dataclass

PROTON_TOP_MOMENTUM = 275.0  # GeV


@dataclass(frozen=True)
class Ion:
    name: str
    A: int
    Z: int
    spin: float
    # Effective nucleon polarizations P_p, P_n building whole-nucleus
    #   g1A = P_p*g1p + P_n*g1n   (naive denominator of the EMC ratio DR)
    # 7Li: P_p=+0.866, P_n=-0.037 (QMC/VMC, Wiringa 1309.3794 via JLab
    #      E12-14-001) -- verified. 3He: Bissey PRC 65:064317 -- verified.
    # 6Li: UNRESOLVED convention (factor 2.4!): Cloet slides use
    #      P_p=P_n=1/3 (per-nucleon-normalized 2-of-6 dilution); cluster
    #      picture gives ~0.81 whole-nucleus (0.87 P_d x 0.93 D-state;
    #      Schellingerhout PRC 48:2714). We keep 1/3 (conservative) until
    #      resolved with I. Cloet -- plans/04_open_questions.md item 6.
    eff_pol_p: float = 0.0
    eff_pol_n: float = 0.0

    @property
    def N(self) -> int:
        return self.A - self.Z

    @property
    def momentum_per_nucleon_max(self) -> float:
        """Rigidity-limited top momentum per nucleon [GeV]."""
        return PROTON_TOP_MOMENTUM * self.Z / self.A


PROTON = Ion("p", 1, 1, 0.5, eff_pol_p=1.0, eff_pol_n=0.0)
DEUTERON = Ion("d", 2, 1, 1.0, eff_pol_p=0.93, eff_pol_n=0.93)  # 1-1.5*w_D
HE3 = Ion("3He", 3, 2, 0.5, eff_pol_p=-0.028, eff_pol_n=0.86)
LI6 = Ion("6Li", 6, 3, 1.0, eff_pol_p=1.0 / 3.0, eff_pol_n=1.0 / 3.0)
LI7 = Ion("7Li", 7, 3, 1.5, eff_pol_p=0.866, eff_pol_n=-0.037)

IONS = {i.name: i for i in (PROTON, DEUTERON, HE3, LI6, LI7)}

ELECTRON_ENERGIES = (5.0, 10.0, 18.0)  # GeV


@dataclass(frozen=True)
class BeamConfig:
    electron_energy: float  # GeV
    ion: Ion
    ion_momentum_per_nucleon: float  # GeV

    @property
    def sqrt_s_per_nucleon(self) -> float:
        """sqrt(s) of the electron-nucleon system [GeV] (massless approx)."""
        return (4.0 * self.electron_energy * self.ion_momentum_per_nucleon) ** 0.5

    def label(self) -> str:
        return (
            f"e({self.electron_energy:g}) x {self.ion.name}"
            f"({self.ion_momentum_per_nucleon:g}/u)"
        )


def default_configs(ion_name: str = "7Li") -> list:
    """Reference energy scan: low/mid/top, mirroring ep 5x41, 10x100, 18x275."""
    ion = IONS[ion_name]
    top = ion.momentum_per_nucleon_max
    return [
        BeamConfig(5.0, ion, round(41.0 * ion.Z / ion.A, 1)),
        BeamConfig(10.0, ion, round(100.0 * ion.Z / ion.A, 1)),
        BeamConfig(18.0, ion, round(top, 1)),
    ]
