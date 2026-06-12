#!/usr/bin/env python3
"""Money plot 1 (plan step 1.3.1): polarized EMC effect on 7Li.

Projected statistical precision on the ratio
  DR(x) = (g1A/F1A) / A1_naive,   A1_naive = [P_p g1p + P_n g1n]/F1A
combined over Q2 bins and the three energy settings, overlaid on the
CBT-like (2x unpolarized EMC) and TMT-like (~1x) scenario curves.
The discrimination significance |DR_CBT - DR_TMT| / dDR is the FOM.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import beams, fom
from polli_fastsim.polarized import (ToyG1, cbt_polarized_emc_ratio,
                                     tmt_polarized_emc_ratio, toy_b1,
                                     toy_delta_gluon)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def delta_dr_per_x(ion_name, lumi, min_events=100):
    """Combined-over-(Q2, energies) statistical error on DR vs x."""
    g1m = ToyG1()
    inv2_tot = None
    x_centers = None
    for cfg in beams.default_configs(ion_name):
        sc = fom.Scenario(lumi_fb_per_nucleon=lumi)
        proj = fom.project_rates(cfg, sc)
        obs = fom.project_observables(cfg, sc, proj, g1m, toy_b1,
                                      toy_delta_gluon)
        x_centers = proj.x[:, 0]
        # A1_naive per bin (no medium modification in the denominator)
        f1 = proj.extras["nf2"].f1a(proj.x, proj.q2) / cfg.ion.A
        g1_naive = g1m.g1_nucleus(cfg.ion, proj.x, proj.q2) / cfg.ion.A
        a1_naive = np.abs(g1_naive / np.maximum(f1, 1e-30))
        err_dr = obs["err_g1_over_f1"] / np.maximum(a1_naive, 1e-12)
        use = proj.accepted & (proj.n_events >= min_events)
        inv2 = np.where(use & (err_dr > 0), 1.0 / err_dr**2, 0.0)
        inv2_tot = inv2.sum(axis=1) if inv2_tot is None \
            else inv2_tot + inv2.sum(axis=1)
    err = np.full(x_centers.shape, np.inf)
    np.divide(1.0, np.sqrt(inv2_tot), out=err, where=inv2_tot > 0)
    return x_centers, err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ion", default="7Li", choices=["6Li", "7Li"])
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    xs = np.logspace(np.log10(0.005), np.log10(0.9), 300)
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7, 7), sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1]})
    ax.plot(xs, cbt_polarized_emc_ratio(xs), "crimson",
            label="CBT-like (2$\\times$ unpol. EMC) [scenario]")
    ax.plot(xs, tmt_polarized_emc_ratio(xs), "navy", ls="--",
            label="TMT-like ($\\approx$1$\\times$) [scenario]")
    ax.axhline(1.0, color="gray", lw=0.5)

    sep = np.abs(cbt_polarized_emc_ratio(xs) - tmt_polarized_emc_ratio(xs))
    for lumi, color, dx in ((10.0, "black", 1.0), (100.0, "seagreen", 1.04)):
        x_c, err = delta_dr_per_x(args.ion, lumi)
        ok = np.isfinite(err) & (x_c > 0.004) & (x_c < 0.9) & (err < 0.5)
        mid = 0.5 * (cbt_polarized_emc_ratio(x_c[ok])
                     + tmt_polarized_emc_ratio(x_c[ok]))
        ax.errorbar(x_c[ok] * dx, mid, yerr=err[ok], fmt=".", ms=4,
                    color=color, lw=1, capsize=0,
                    label=f"proj. $\\delta\\Delta R$, {lumi:g} fb$^{{-1}}$/u")
        sep_c = np.abs(cbt_polarized_emc_ratio(x_c[ok])
                       - tmt_polarized_emc_ratio(x_c[ok]))
        ax2.plot(x_c[ok], sep_c / err[ok], ".-", color=color, lw=1, ms=4,
                 label=f"{lumi:g} fb$^{{-1}}$/u")

    ax.set_xscale("log")
    ax.set_ylim(0.5, 1.3)
    ax.set_ylabel(r"$\Delta R(x) = g_1^A/(P_p g_1^p + P_n g_1^n)$")
    ax.set_title(
        f"Polarized EMC effect, {args.ion}: scenario discrimination\n"
        "(3 energy settings combined; $P_e$=0.7, $P_z$=0.7; "
        "stat. only; TOY inputs)", fontsize=10)
    ax.legend(fontsize=8, loc="lower left")
    ax2.axhline(5, color="gray", ls=":", lw=1)
    ax2.set_yscale("log")
    ax2.set_xlabel(r"$x$")
    ax2.set_ylabel(r"$|\Delta R_{CBT}-\Delta R_{TMT}|/\delta$")
    ax2.legend(fontsize=8)
    fig.tight_layout()
    path = outdir / f"money_polemc_{args.ion}.png"
    fig.savefig(path, dpi=150)
    print(f"wrote {path}")
    x_c, err10 = delta_dr_per_x(args.ion, 10.0)
    for xq in (0.1, 0.3, 0.5, 0.7):
        i = np.argmin(np.abs(x_c - xq))
        print(f"  x={x_c[i]:.2f}: dDR(10/fb)={err10[i]:.4f}  "
              f"CBT-TMT sep={abs(cbt_polarized_emc_ratio(x_c[i]) - tmt_polarized_emc_ratio(x_c[i])):.3f}")


if __name__ == "__main__":
    main()
