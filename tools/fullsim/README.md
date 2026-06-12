# Full-simulation runbook (Phase 2)

Local ePIC chain status and recipes. Versions matter — record them in
every output directory.

## Containers

| container | state |
|---|---|
| `~/Projects/eic/local/lib/jug_xl-nightly.sif` | legacy (frozen ~Sep 2024); detector configs `/opt/detector/epic-24.0{5,6,7}.0`; works — used for the first gun scan |
| `~/Projects/eic/local/lib/eic_xl-nightly.sif` | **current** (installed 2026-06-12): EICrecon v1.38.0, detector configs epic-26.03.0…26.06.0 + `epic-main`; source `/opt/detector/epic-main/bin/thisepic.sh` |

`./eic-shell` now starts the new container. Re-run the gun scan in it
(epic-26.06 far-forward geometry) before quoting acceptance numbers.

**Known issue (found 2026-06-12):** the `eic_xl-nightly` pyHepMC3
`rootIO.ReaderRootTree` **segfaults** on files the legacy container reads
fine — use `jug_xl-nightly.sif` for all HepMC3 tree.root *reading*
(dump_spectators.py); npsim/eicrecon in the new container are unaffected.
Consider reporting upstream (eic-shell / pyHepMC3 bindings).

Inside the legacy container, source the geometry with:
```bash
source /opt/detector/epic-24.07.0/setup.sh   # sets $DETECTOR_PATH
```
(current containers use `/opt/detector/epic-main/bin/thisepic.sh`).

## Far-forward gun scan (plan 2.2.1)

```bash
SIF=~/Projects/eic/local/lib/jug_xl-nightly.sif
singularity exec $SIF bash ff_gun_scan.sh ~/Projects/eic/data/ff_gun 60
singularity exec $SIF python3 ff_gun_hits.py ~/Projects/eic/data/ff_gun
```

The scan shoots fragments along the ion axis (−25 mrad in x) with total
momenta chosen so the rigidity ratio R relative to the default
craterlake fields (18×275 ep optics) reproduces the Li-fragment cases:

| config | particle | p [GeV] | R | expectation (plans/03 §2.2) |
|---|---|---|---|---|
| alpha_R0857 | α | 471.4 | 0.857 | Roman Pots (⁷Li α-tag) |
| alpha_R100 | α | 550.0 | 1.000 | lost in beam pipe (⁶Li α, pT≈0) |
| alpha_R100pt3 | α | 550.0 | 1.000 (+0.55 mrad ≙ pT=0.3) | RP pT-tail |
| triton_R1286 | t | 353.6 | 1.286 | no coverage |
| proton_R050 | p | 137.5 | 0.500 | OMD |
| deuteron_R0857 | d | 235.7 | 0.857 | Roman Pots |
| neutron_ZDC | n | 117.9 | — | ZDC |

**Results (2026-06-12, epic-main geometry, 18×275 fields, 12 ev/config —
fraction of events with ≥1 hit):**

| config | RP | OMD | B0 | ZDC | verdict vs prediction |
|---|---|---|---|---|---|
| alpha_R0857 | **1.00** | 0 | 0 | 0.17 | ✓ RP (⁷Li α-tag works) |
| alpha_R100 | 0.08 | 0 | 0 | 0 | ✓ invisible (⁶Li α beam-blind) |
| alpha_R100pt3 | **1.00** | 0 | 0 | 0 | ✓ pT-tail recovers it |
| deuteron_R0857 | **1.00** | 0 | 0 | 0.25 | ✓ RP |
| proton_R050 | 0.25 | **1.00** | 0 | 0.42 | ✓ OMD (RP/ZDC = pass-through/splash) |
| triton_R1286 | **1.00** | 0 | 0 | **1.00** | ✗ **surprise: taggable!** |
| neutron_ZDC | 0 | 0 | 0 | **1.00** | ✓ ZDC |

Findings beyond the table:
- Current epic-main RP stations sit at **z ≈ 32.5 / 34.3 m** (papers
  quote 26/28 m — geometry has moved; update plans/03 numbers when the
  preTDR layout is confirmed).
- Hit x-positions: α(R=0.857) crosses ~40 mm on the dispersion side of
  the beam orbit; the **over-rigid triton (R=1.286) crosses ~36 mm on
  the opposite side of the *same* RP planes and then deposits in the
  ZDC** (it bends less than the beam, staying near the zero-degree
  line). The "no IP6 triton coverage" assumption may be wrong — the
  ⁷Li α+t double-tag could work with RP-inner-side tracking + ZDC
  energy. Needs: 10σ beam-exclusion modeling, divergence/vertex spread,
  reconstruction-level confirmation.
- Gotcha encoded in the script: plain `epic_craterlake.xml` loads the
  **5×41 beamline fields** — with the 275-optics momenta everything
  flies straight to the ZDC. Use `epic_craterlake_18x275.xml`.
- `epic_craterlake_18x110_Au.xml` (Z/A ≈ 0.40) is the closest existing
  optics proxy for a ⁷Li beam (Z/A = 0.429).

Caveats: single-particle, fixed vertex/direction (no beam envelope —
real pots retract to 10σ), hit-level only. This validates *routing
logic*, not final acceptances (plans/03 step 2.2 proper).

## e+d control inputs (plan 1.5.3)

Official BeAGLE samples streamed via xrootd (no FLUKA needed):
```bash
singularity exec $SIF python3 ../analysis/dump_spectators.py \
  'root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-1.3/eH2/en/10x130/q2_1to1000/..._run001.hepmc3.tree.root' \
  out.csv --nevents 50000
python3 ../analysis/ed_control_analysis.py out.csv
```
