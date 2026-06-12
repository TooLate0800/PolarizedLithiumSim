#!/usr/bin/env python3
"""Count far-forward hits per subsystem for the gun-scan outputs.
Run INSIDE eic-shell (needs uproot):
  singularity exec $SIF python3 ff_gun_hits.py <outdir-from-ff_gun_scan>
"""

import glob
import re
import sys

import uproot

PATTERNS = {
    "RomanPots": re.compile(r"ForwardRomanPot.*Hits$"),
    "OMD": re.compile(r"ForwardOffM.*Hits$"),
    "B0": re.compile(r"B0Tracker.*Hits$|B0ECal.*Hits$"),
    "ZDC": re.compile(r"(Hcal|Ecal)FarForwardZDC.*Hits$|ZDC.*Hits$"),
}


def main():
    outdir = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ff_gun"
    files = sorted(glob.glob(f"{outdir}/gun_*.edm4hep.root"))
    if not files:
        sys.exit(f"no gun_*.edm4hep.root in {outdir}")
    print(f"{'config':16s} {'nev':>4s} " +
          " ".join(f"{k:>10s}" for k in PATTERNS) + "   (fraction of events with >=1 hit)")
    for f in files:
        tag = f.split("gun_")[-1].replace(".edm4hep.root", "")
        try:
            t = uproot.open(f)["events"]
        except Exception as e:
            print(f"{tag:16s} OPEN-FAIL {e}")
            continue
        nev = t.num_entries
        cols = {}
        for label, pat in PATTERNS.items():
            frac = 0.0
            for b in t.keys():
                base = b.split("/")[0].split(".")[0]
                if pat.match(base) and b.endswith(".cellID"):
                    counts = t[b].array(library="np")
                    frac = max(frac,
                               float(sum(1 for c in counts if len(c) > 0))
                               / max(nev, 1))
            cols[label] = frac
        print(f"{tag:16s} {nev:4d} " +
              " ".join(f"{cols[k]:10.2f}" for k in PATTERNS))
    # list available FF-ish collections once, for reference
    t = uproot.open(files[0])["events"]
    ff = sorted({b.split(".")[0] for b in t.keys()
                 if re.search(r"RomanPot|OffM|B0|ZDC", b)})
    print("\ncollections present:", ", ".join(ff[:12]))


if __name__ == "__main__":
    main()
