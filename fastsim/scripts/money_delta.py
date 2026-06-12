#!/usr/bin/env python3
"""Money plot 3 (plan step 1.3.3): gluonometry discovery reach.

5-sigma luminosity for the cos(2phi) double-helicity-flip amplitude as a
function of the Delta/F1 scale (scenario range 1e-3 .. 1e-2 spans the
Sather-Schmidt bag estimate). Assumes transversely tensor-polarized 6Li
running with unpolarized electrons; significance combines all accepted
(x, Q2) bins:  sig^2 = sum_bins A_bin^2 * P_zz^2 * N_bin / 2.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import a_cos2phi
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def lumi_5sigma(cfg, scale, pzz, min_events=10):
    """fb^-1/nucleon for 5 sigma, combining all accepted bins."""
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz)
    proj = fom.project_rates(cfg, sc)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    use = proj.accepted & (proj.n_events >= min_events)
    sig2_per_fb = np.where(use,
                           amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()
    return 25.0 / max(sig2_per_fb, 1e-30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="6Li", choices=["6Li", "7Li"])
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    scales = np.logspace(-3.3, -1.7, 15)
    fig, ax = plt.subplots(figsize=(7, 5))
    for cfg, color in zip(beams.default_configs(args.ion),
                          ("crimson", "seagreen", "navy")):
        for pzz, ls in ((0.60, "-"), (0.80, "--")):
            reach = [lumi_5sigma(cfg, s, pzz) for s in scales]
            ax.plot(scales, reach, ls, color=color, lw=1.5,
                    label=f"{cfg.label()}, $P_{{zz}}$={pzz:g}")
    ax.axhspan(1, 100, color="gold", alpha=0.12,
               label="1-100 fb$^{-1}$/u (plausible program)")
    ax.axvline(1e-3, color="gray", ls=":", lw=1)
    ax.text(1.05e-3, 0.93, "Sather-Schmidt\n$O(10^{-3})$", fontsize=7,
            transform=ax.get_xaxis_transform(), va="top")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta/F_1$ scale (peak of scenario shape)")
    ax.set_ylabel(r"$L_{5\sigma}$ [fb$^{-1}$/nucleon]")
    ax.set_title(f"Nuclear gluonometry reach, transversely polarized "
                 f"{args.ion}\n(cos 2$\\phi$ amplitude, all bins combined; "
                 "TOY inputs)", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    path = outdir / f"money_delta_{args.ion}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    for cfg in beams.default_configs(args.ion):
        print(f"  {cfg.label():26s} L_5sig(Delta/F1=1e-3, Pzz=0.8) = "
              f"{lumi_5sigma(cfg, 1e-3, 0.8):9.1f} fb^-1/u ; "
              f"(1e-2) = {lumi_5sigma(cfg, 1e-2, 0.8):7.3f}")


if __name__ == "__main__":
    main()
