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
from polli_fastsim.inputs import get_backends
from polli_fastsim.polarized import toy_delta_gluon
from polli_fastsim.structure import NuclearF2

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def sig2_per_fb_at(cfg, scale, pzz, base=None, min_events=10):
    """Significance^2 per fb^-1/nucleon at a reference Delta/F1 scale.
    The toy Delta shape is linear in `scale`, so reach curves follow
    L_5sig(s) = 25 / (sig2 * (s/scale)^2) analytically."""
    sc = fom.Scenario(lumi_fb_per_nucleon=1.0, pol_ion_tensor=pzz)
    nf2_in = NuclearF2(cfg.ion, base=base) if base is not None else None
    proj = fom.project_rates(cfg, sc, nuclear_f2=nf2_in)
    nf2 = proj.extras["nf2"]
    f1 = nf2.f1a(proj.x, proj.q2) / cfg.ion.A
    f2 = nf2.f2a(proj.x, proj.q2) / cfg.ion.A
    y = proj.extras["y"]
    delta = toy_delta_gluon(proj.x, proj.q2, f1, scale=scale)
    amp = a_cos2phi(delta, f1, f2, proj.x, y)
    use = proj.accepted & (proj.n_events >= min_events)
    return np.where(use, amp**2 * pzz**2 * proj.n_events / 2.0, 0.0).sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="6Li", choices=["6Li", "7Li"])
    ap.add_argument("--pdf", default="toy", choices=["toy", "grid"])
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    backends = get_backends(args.pdf)
    scales = np.logspace(-3.3, -1.7, 15)
    s0 = 1e-3
    fig, ax = plt.subplots(figsize=(7, 5))
    sig2_ref = {}
    for cfg, color in zip(beams.default_configs(args.ion),
                          ("crimson", "seagreen", "navy")):
        for pzz, ls in ((0.60, "-"), (0.80, "--")):
            sig2 = sig2_per_fb_at(cfg, s0, pzz, base=backends["base"])
            sig2_ref[(cfg.label(), pzz)] = sig2
            reach = 25.0 / np.maximum(sig2 * (scales / s0) ** 2, 1e-30)
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
                 f"{backends['tag'].upper()} inputs)", fontsize=10)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    path = outdir / f"money_delta_{args.ion}_{backends['tag']}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    for cfg in beams.default_configs(args.ion):
        sig2 = sig2_ref[(cfg.label(), 0.8)]
        print(f"  {cfg.label():26s} L_5sig(Delta/F1=1e-3, Pzz=0.8) = "
              f"{25.0/sig2:9.1f} fb^-1/u ; (1e-2) = {25.0/sig2/100:7.3f}")


if __name__ == "__main__":
    main()
