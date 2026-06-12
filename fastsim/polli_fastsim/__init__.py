"""Fast-simulation toolkit for the polarized 6,7Li @ EIC feasibility study.

Phase-1 starter package: kinematic phase space, rate estimates, and
statistical figure-of-merit projections for the three flagship observables
of the polarized-lithium program (docs/ecrp_2026_proposal.pdf):

  1. Polarized EMC effect via g1 of 7Li (and 6Li)        -> A_parallel
  2. Tensor structure function b1 of spin-1 6Li           -> A_zz
  3. Gluon transversity Delta(x,Q2) "nuclear gluonometry" -> A_cos(2phi)

All structure-function inputs are intentionally swappable: the default
"toy" parameterizations are good to tens of percent for rate/coverage
studies and MUST be replaced by LHAPDF/JAM/DSSV grids (and theory curves
for b1/Delta) before quoting figures of merit in any document.
"""

from . import beams, kinematics, structure, polarized, asymmetries, fom

__all__ = ["beams", "kinematics", "structure", "polarized", "asymmetries", "fom"]
__version__ = "0.1.0"
