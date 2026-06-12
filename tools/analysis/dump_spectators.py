#!/usr/bin/env python3
"""Stream a (possibly remote) HepMC3 tree.root file and dump spectator
candidates + beams + scattered electron to CSV. Runs INSIDE eic-shell:

  singularity exec $SIF python3 dump_spectators.py \
      root://dtn-eic.jlab.org//volatile/.../file.hepmc3.tree.root \
      out.csv --nevents 50000

Rows: ievt, kind, pdg, px, py, pz, e   (kind: B=beam-ion, E=scattered-e,
S=baryon/nucleus candidate with |eta|>4 or pz>0.1*beam)
"""

import argparse
import csv
import sys

from pyHepMC3 import HepMC3
from pyHepMC3.rootIO import HepMC3 as rootIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--nevents", type=int, default=50000)
    args = ap.parse_args()

    reader = rootIO.ReaderRootTree(args.input)
    if reader.failed():
        sys.exit(f"cannot open {args.input}")
    out = csv.writer(open(args.output, "w"))
    out.writerow(["ievt", "kind", "pdg", "px", "py", "pz", "e"])

    evt = HepMC3.GenEvent()
    n = 0
    while n < args.nevents and not reader.failed():
        if not reader.read_event(evt) or reader.failed():
            break
        beam_pz = 1.0
        rows = []
        scat_e = None
        for p in evt.particles():
            st, pid = p.status(), p.pid()
            mom = p.momentum()
            if st == 4 and pid != 11:
                rows.append((n, "B", pid, mom.px(), mom.py(), mom.pz(),
                             mom.e()))
                beam_pz = max(beam_pz, mom.pz())
            elif st == 1:
                if pid == 11:
                    if scat_e is None or mom.e() > scat_e[6]:
                        scat_e = (n, "E", pid, mom.px(), mom.py(), mom.pz(),
                                  mom.e())
                elif pid in (2212, 2112) or pid > 1000000000:
                    rows.append((n, "S", pid, mom.px(), mom.py(), mom.pz(),
                                 mom.e()))
        # keep only forward (ion-side) baryons: pz > 10% of beam pz
        for r in rows:
            if r[1] == "B" or r[5] > 0.1 * beam_pz:
                out.writerow(r)
        if scat_e is not None:
            out.writerow(scat_e)
        n += 1
        if n % 10000 == 0:
            print(f"  {n} events", flush=True)
    reader.close()
    print(f"dumped {n} events to {args.output}")


if __name__ == "__main__":
    main()
