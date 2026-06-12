#!/usr/bin/env python3
"""First quantitative spectator-tagging acceptance estimates for e+6,7Li.

Cluster-model spectators (polli_fastsim.spectator) folded with the
verified far-forward windows (polli_fastsim.farforward). The headline
physics question (plans/02 step 1.5.4): does the soft alpha from 6Li
survive the Roman-Pot near-beam pT cutoff, and how cleanly does the
off-rigidity alpha from 7Li land in the Roman Pots?

Outputs: acceptance tables (text) for both optics and a beta-scan of the
wave-function tail, plus pT / R / theta distributions per channel.
"""

import argparse
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ENERGIES = {"6Li": 137.5, "7Li": 117.9}  # GeV/u top rigidity settings


def channel_plots(channel, kin, outdir, tag):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6))
    axes[0].hist(kin["pT"], bins=120, range=(0, 0.6), histtype="step",
                 density=True)
    for cut, ls in ((0.20, "--"), (0.45, ":")):
        axes[0].axvline(cut, color="crimson", ls=ls, lw=1,
                        label=f"RP pT cut {cut:g} GeV")
    axes[0].set_xlabel(r"$p_T$ [GeV]")
    axes[0].legend(fontsize=7)
    axes[1].hist(kin["R"], bins=120, histtype="step", density=True)
    for lo, hi, c in (( *ff.RP_R_WINDOW, "seagreen"),
                      (*ff.OMD_R_WINDOW, "orange")):
        axes[1].axvspan(lo, hi, alpha=0.15, color=c)
    axes[1].set_xlabel(r"rigidity ratio $R$")
    axes[2].hist(kin["theta"] * 1e3, bins=120, histtype="step", density=True)
    axes[2].axvline(5.0, color="gray", ls="--", lw=1)
    axes[2].set_xlabel(r"$\theta$ [mrad]")
    fig.suptitle(channel.name, fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f"spectator_{tag}.png", dpi=140)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nevents", type=int, default=400_000)
    ap.add_argument("--outdir", default="out")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    lines = ["# Spectator-tagging acceptance (cluster model, far-forward "
             "windows of plans/03 §2.2)",
             "# beta = wave-function short-range scale [GeV] (tail "
             "uncertainty scan)", ""]
    rng = np.random.default_rng(7)

    for channel in sp.CHANNELS:
        isotope = "6Li" if channel.beam_A == 6 else "7Li"
        p_u = ENERGIES[isotope]
        lines.append(f"== {channel.name}  (beam {p_u:g} GeV/u, "
                     f"kappa = {channel.kappa*1e3:.1f} MeV, "
                     f"L = {channel.l_wave}) ==")
        for beta in (0.20, 0.30, 0.40):
            kin = sp.spectator_lab_kinematics(channel, p_u, args.nevents,
                                              beta=beta, rng=rng)
            for optics in (ff.HIGH_ACCEPTANCE, ff.HIGH_DIVERGENCE):
                acc = ff.acceptance_summary(kin["R"], kin["theta"],
                                            kin["pT"], optics)
                tagged = 1.0 - acc["lost"]
                parts = "  ".join(f"{k}:{v:6.3f}" for k, v in acc.items()
                                  if v > 5e-4 and k != "lost")
                lines.append(f"  beta={beta:.2f}  {optics.name:16s} "
                             f"tagged={tagged:6.3f}   {parts}")
            if beta == 0.30:
                tag = (channel.spectator + "_" + isotope).replace("/", "")
                channel_plots(channel, kin, outdir, tag)
        lines.append("")

    text = "\n".join(lines)
    (outdir / "tagging_acceptance.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
