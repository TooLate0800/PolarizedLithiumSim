#!/usr/bin/env python3
"""e+d control study (plan step 1.5.3): official BeAGLE sample vs the
cluster-spectator model + far-forward routing.

Input: CSV from dump_spectators.py on a BeAGLE eH2 'en' sample (struck
neutron, spectator proton). Beam row gives the event-by-event ion axis
(afterburned: -25 mrad crossing + divergence); spectator kinematics are
computed RELATIVE to that axis, which is also what the far-forward optics
see.

Outputs: pT and xL spectra overlaid with the deuteron Hulthen model,
far-forward routing fractions (BeAGLE vs model), summary text.
"""

import argparse
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "fastsim"))

from polli_fastsim import farforward as ff
from polli_fastsim import spectator as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load(csv_path):
    raw = np.genfromtxt(csv_path, delimiter=",", names=True, dtype=None,
                        encoding="utf-8")
    kinds = raw["kind"].astype(str)
    return raw, kinds


def spectator_protons(raw, kinds, pdg=2212, beam_a_over_z=2.0):
    """Per event: beam axis from 'B' row; spectator = baryon (given pdg)
    with largest pz; returns (pT_rel, xL, R, theta) arrays."""
    beams = raw[kinds == "B"]
    prot = raw[(kinds == "S") & (raw["pdg"] == pdg)]
    beam_by_evt = {int(b["ievt"]): b for b in beams}
    best = {}
    for p in prot:
        i = int(p["ievt"])
        if i not in beam_by_evt:
            continue
        if i not in best or p["pz"] > best[i]["pz"]:
            best[i] = p
    pt, xl, rr, th = [], [], [], []
    for i, p in best.items():
        b = beam_by_evt[i]
        bv = np.array([b["px"], b["py"], b["pz"]])
        pn = np.linalg.norm(bv)          # per-nucleon beam momentum
        u = bv / pn
        pv = np.array([p["px"], p["py"], p["pz"]])
        p_par = float(pv @ u)
        p_tot = float(np.linalg.norm(pv))
        p_perp = float(np.sqrt(max(p_tot**2 - p_par**2, 0.0)))
        if p_par < 0.6 * pn:             # spectator selection: xL > 0.6
            continue
        pt.append(p_perp)
        xl.append(p_par / pn)
        # beam record is per-nucleon: beam rigidity = pn * A_beam/Z_beam
        rr.append(p_tot / (pn * beam_a_over_z))
        th.append(np.arctan2(p_perp, p_par))
    return map(np.asarray, (pt, xl, rr, th))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--outdir", default="fastsim/out")
    ap.add_argument("--beta", type=float, default=0.26,
                    help="Hulthen short-range scale for the model overlay")
    ap.add_argument("--spectator", default="p", choices=["p", "n"])
    ap.add_argument("--beam", default="d", choices=["d", "He3"])
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    pdg = 2212 if args.spectator == "p" else 2112
    if args.beam == "He3":
        channel = sp.HE3_P_TAG
    else:
        channel = (sp.DEUTERON_P_TAG if args.spectator == "p"
                   else sp.DEUTERON_N_TAG)
    raw, kinds = load(args.csv)
    aoz = 1.5 if args.beam == "He3" else 2.0
    pt, xl, rr, th = spectator_protons(raw, kinds, pdg=pdg,
                                       beam_a_over_z=aoz)
    n_evt = len(set(raw["ievt"].astype(int)))
    lines = [f"# e+d control: {args.csv}",
             f"events: {n_evt}, spectator-{args.spectator} candidates: "
             f"{len(pt)} ({len(pt)/max(n_evt,1):.1%} of events)"]

    pn_model = 166.0 if args.beam == "He3" else 130.0
    kin = sp.spectator_lab_kinematics(channel, pn_model, 300000,
                                      beta=args.beta)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    axes[0].hist(pt, bins=80, range=(0, 0.8), density=True, histtype="step",
                 color="black",
                 label=f"BeAGLE e+{args.beam} (official sample)")
    axes[0].hist(kin["pT"], bins=80, range=(0, 0.8), density=True,
                 histtype="step", color="crimson",
                 label=f"Hulthen cluster model (beta={args.beta:g})")
    axes[0].set_xlabel(r"spectator-p $p_T$ w.r.t. ion axis [GeV]")
    axes[0].set_yscale("log")
    axes[0].legend(fontsize=7)
    axes[1].hist(xl, bins=80, range=(0.7, 1.3), density=True,
                 histtype="step", color="black")
    axes[1].hist(kin["xL"], bins=80, range=(0.7, 1.3), density=True,
                 histtype="step", color="crimson")
    axes[1].set_xlabel(r"$x_L = p_\parallel / p_{beam/u}$")
    if args.spectator == "p":
        rlo, rhi = (0.45, 0.85) if args.beam == "He3" else (0.35, 0.65)
        axes[2].hist(rr, bins=80, range=(rlo, rhi), density=True,
                     histtype="step", color="black")
        axes[2].hist(kin["R"], bins=80, range=(rlo, rhi), density=True,
                     histtype="step", color="crimson")
        for lo, hi, c in ((*ff.OMD_R_WINDOW, "orange"),
                          (*ff.RP_R_WINDOW, "g")):
            axes[2].axvspan(lo, hi, alpha=0.12, color=c)
        axes[2].set_xlabel(f"rigidity ratio R ({args.beam} beam)")
    else:
        axes[2].hist(th * 1e3, bins=80, range=(0, 6), density=True,
                     histtype="step", color="black")
        axes[2].hist(kin["theta"] * 1e3, bins=80, range=(0, 6),
                     density=True, histtype="step", color="crimson")
        axes[2].axvline(ff.THETA_ZDC_MAX * 1e3, color="gray", ls="--")
        axes[2].set_xlabel(r"$\theta$ w.r.t. ion axis [mrad] (ZDC window)")
    fig.suptitle(f"e+{args.beam} spectator-{args.spectator} control: "
                 "official BeAGLE vs cluster model (afterburned)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f"ed_control_spectra_{args.beam}_{args.spectator}.png",
                dpi=140)

    if args.spectator == "p":
        pairs = (("BeAGLE", (rr, th, pt)),
                 ("model", (kin["R"], kin["theta"], kin["pT"])))
        for name, (R, T, P) in pairs:
            for optics in (ff.HIGH_ACCEPTANCE, ff.HIGH_DIVERGENCE):
                acc = ff.acceptance_summary(R, T, P, optics)
                parts = "  ".join(f"{k}:{v:6.3f}" for k, v in acc.items()
                                  if v > 5e-4)
                lines.append(f"{name:7s} {optics.name:16s} {parts}")
    else:
        for name, T in (("BeAGLE", th), ("model", kin["theta"])):
            acc = ff.neutral_summary(T)
            lines.append(f"{name:7s} ZDC:{acc['ZDC']:6.3f}  "
                         f"lost:{acc['lost']:6.3f}")
    # quantitative tail comparison
    for cut in (0.1, 0.2, 0.3, 0.45):
        fb = float(np.mean(pt > cut))
        fm = float(np.mean(kin["pT"] > cut))
        lines.append(f"P(pT > {cut:4.2f}) BeAGLE={fb:6.4f}  model={fm:6.4f}"
                     f"  ratio={fb/max(fm,1e-9):5.2f}")
    text = "\n".join(lines)
    (outdir / f"ed_control_summary_{args.beam}_{args.spectator}.txt"
     ).write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
