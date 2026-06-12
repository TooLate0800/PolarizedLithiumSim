#!/usr/bin/env python3
"""Money plot 2 (plan step 1.3.2): tensor asymmetry A_zz / b1 on 6Li.

Projected per-x statistical precision on A_zz (longitudinal tensor,
equal-thirds spin states) vs the two scenario curves: 'HERMES-like'
(large low-x b1) and 'standard convolution' (10x smaller), both scaled
to the embedded deuteron in 6Li (x P_d = 0.87). Shown for two tensor
polarizations -- this plot sets the P_zz requirement for the source.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.asymmetries import azz as azz_formula
from polli_fastsim.polarized import (ToyG1, b1_convolution,
                                     b1_li6_from_deuteron, toy_b1,
                                     toy_delta_gluon)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lumi", type=float, default=10.0)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    g1m = ToyG1()
    fig, ax = plt.subplots(figsize=(7, 5))

    # scenario Azz(x) curves at a representative Q2 slice (combined plot
    # uses the bin-by-bin values; curves drawn at Q2 = 4 GeV2)
    xs = np.logspace(np.log10(0.003), np.log10(0.7), 200)
    q2s = np.full_like(xs, 4.0)
    cfg0 = beams.default_configs("6Li")[0]
    s0 = cfg0.sqrt_s_per_nucleon**2
    y0 = np.clip(q2s / (s0 * xs), 1e-4, 1.0)
    from polli_fastsim.structure import NuclearF2
    nf2 = NuclearF2(cfg0.ion)
    f1 = nf2.f1a(xs, q2s) / cfg0.ion.A
    f2 = nf2.f2a(xs, q2s) / cfg0.ion.A
    for b1f, color, label in (
            (toy_b1, "crimson", "HERMES-like $b_1$ [scenario]"),
            (b1_convolution, "navy", "convolution $b_1$ [scenario]")):
        b1 = b1_li6_from_deuteron(b1f(xs, q2s, f1))
        ax.plot(xs, np.abs(azz_formula(b1, f1, f2, xs, y0)), color=color,
                label=label)

    # projected per-x errors, combined over Q2 and the 3 energy settings
    for pzz, fmt in ((0.60, "ko"), (0.80, "g^")):
        inv2_tot, x_c = None, None
        for cfg in beams.default_configs("6Li"):
            sc = fom.Scenario(lumi_fb_per_nucleon=args.lumi,
                              pol_ion_tensor=pzz)
            proj = fom.project_rates(cfg, sc)
            obs = fom.project_observables(cfg, sc, proj, g1m, toy_b1,
                                          toy_delta_gluon)
            x_c = proj.x[:, 0]
            use = proj.accepted & (proj.n_events >= 100)
            inv2 = np.where(use, 1.0 / obs["err_azz"]**2, 0.0).sum(axis=1)
            inv2_tot = inv2 if inv2_tot is None else inv2_tot + inv2
        err = np.full(x_c.shape, np.inf)
        np.divide(1.0, np.sqrt(inv2_tot), out=err, where=inv2_tot > 0)
        ok = np.isfinite(err) & (x_c > 0.003) & (x_c < 0.7)
        ax.plot(x_c[ok], err[ok], fmt, ms=4, ls="-", lw=0.8,
                label=f"$\\delta A_{{zz}}$/x-bin, $P_{{zz}}$={pzz:g}, "
                      f"{args.lumi:g} fb$^{{-1}}$/u")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$|A_{zz}|$  and  $\delta A_{zz}$")
    ax.set_ylim(1e-5, 0.3)
    ax.set_title("Tensor asymmetry on $^6$Li (embedded deuteron, "
                 "$P_d$=0.87)\n3 energies combined; stat. only; TOY inputs",
                 fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / "money_b1_6Li.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
