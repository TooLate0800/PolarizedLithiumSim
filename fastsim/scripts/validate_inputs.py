#!/usr/bin/env python3
"""Step 1.2 validation: toy structure functions vs LHAPDF grids.

Compares ToyF2 against CT18NLO (LO F2 formula) and ToyG1's A1p against
NNPDFpol11_100, at representative (x, Q2) points. Documents the accuracy
of the toy inputs and certifies the grid backends.
"""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim.polarized import PartonG1, ToyG1
from polli_fastsim.structure import PartonF2, ToyF2

POINTS = [(1e-3, 10.0), (1e-2, 10.0), (0.1, 10.0), (0.3, 10.0),
          (0.5, 10.0), (0.7, 25.0), (0.01, 100.0), (0.1, 1000.0)]


def main():
    toy_f2, grid_f2 = ToyF2(), PartonF2()
    toy_g1, grid_g1 = ToyG1(), PartonG1()
    print(f"{'x':>8} {'Q2':>7} | {'F2p toy':>8} {'F2p CT18':>9} {'ratio':>6}"
          f" | {'A1p toy':>8} {'A1p NNpol':>9}")
    worst = 1.0
    for x, q2 in POINTS:
        f2t = float(toy_f2.f2p(x, q2))
        f2g = float(grid_f2.f2p(x, q2))
        r = f2t / f2g if f2g > 0 else np.nan
        worst = max(worst, max(r, 1 / r) if r > 0 else worst)
        a1t = float(toy_g1.a1p(x))
        f1g = f2g / (2 * x * 1.18)  # crude R for the comparison only
        a1g = float(grid_g1.g1p(x, q2)) / max(f1g, 1e-30)
        print(f"{x:8.3g} {q2:7.1f} | {f2t:8.3f} {f2g:9.3f} {r:6.2f}"
              f" | {a1t:8.3f} {a1g:9.3f}")
    print(f"\nworst toy/grid F2 ratio: {worst:.2f} "
          "(toy is for coverage maps only; FOM scripts can switch to "
          "PartonF2/PartonG1 with one line)")


if __name__ == "__main__":
    main()
