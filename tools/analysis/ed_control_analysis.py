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


def spectator_protons(raw, kinds):
    """Per event: beam axis from 'B' row; spectator = proton with largest
    pz; returns (pT_rel, xL, R, theta) arrays."""
    beams = raw[kinds == "B"]
    prot = raw[(kinds == "S") & (raw["pdg"] == 2212)]
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
        # deuteron beam rigidity = 2*pn/Z(=1); proton rigidity = p_tot/1
        rr.append(p_tot / (2.0 * pn))
        th.append(np.arctan2(p_perp, p_par))
    return map(np.asarray, (pt, xl, rr, th))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--outdir", default="fastsim/out")
    ap.add_argument("--beta", type=float, default=0.26,
                    help="Hulthen short-range scale for the model overlay")
    args = ap.parse_args()
    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    raw, kinds = load(args.csv)
    pt, xl, rr, th = spectator_protons(raw, kinds)
    n_evt = len(set(raw["ievt"].astype(int)))
    lines = [f"# e+d control: {args.csv}",
             f"events: {n_evt}, spectator-proton candidates: {len(pt)} "
             f"({len(pt)/max(n_evt,1):.1%} of events)"]

    # model: deuteron Hulthen, proton spectator, beam 130 GeV/u
    kin = sp.spectator_lab_kinematics(sp.DEUTERON_P_TAG, 130.0, 300000,
                                      beta=args.beta)

    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    axes[0].hist(pt, bins=80, range=(0, 0.8), density=True, histtype="step",
                 color="black", label="BeAGLE eH2 'en' (official)")
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
    axes[2].hist(rr, bins=80, range=(0.35, 0.65), density=True,
                 histtype="step", color="black")
    axes[2].hist(kin["R"], bins=80, range=(0.35, 0.65), density=True,
                 histtype="step", color="crimson")
    for lo, hi, c in ((*ff.OMD_R_WINDOW, "orange"), (*ff.RP_R_WINDOW, "g")):
        axes[2].axvspan(lo, hi, alpha=0.12, color=c)
    axes[2].set_xlabel("rigidity ratio R (d beam)")
    fig.suptitle("e+d spectator-proton control: official BeAGLE vs cluster "
                 "model (10x130, afterburned)", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "ed_control_spectra.png", dpi=140)

    for name, (R, T, P) in (("BeAGLE", (rr, th, pt)),
                            ("model", (kin["R"], kin["theta"], kin["pT"]))):
        for optics in (ff.HIGH_ACCEPTANCE, ff.HIGH_DIVERGENCE):
            acc = ff.acceptance_summary(R, T, P, optics)
            parts = "  ".join(f"{k}:{v:6.3f}" for k, v in acc.items()
                              if v > 5e-4)
            lines.append(f"{name:7s} {optics.name:16s} {parts}")
    # quantitative tail comparison
    for cut in (0.1, 0.2, 0.3, 0.45):
        fb = float(np.mean(pt > cut))
        fm = float(np.mean(kin["pT"] > cut))
        lines.append(f"P(pT > {cut:4.2f}) BeAGLE={fb:6.4f}  model={fm:6.4f}"
                     f"  ratio={fb/max(fm,1e-9):5.2f}")
    text = "\n".join(lines)
    (outdir / "ed_control_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
