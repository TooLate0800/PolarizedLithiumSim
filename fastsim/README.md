# polli_fastsim — Phase-1 fast simulation starter

Fast (analytic, binned) simulation for the polarized ⁶Li/⁷Li @ EIC
feasibility study. No event generation: rates per (x,Q²) bin from a
structure-function model × luminosity, statistical figure-of-merit from
asymmetry estimators. Companion to `../plans/02_phase1_event_generation.md`.

## Quick start

```bash
cd fastsim
python3 -m pytest tests/ -q          # smoke tests
python3 scripts/phase_space_map.py --ion 7Li --lumi 10 --outdir out
python3 scripts/phase_space_map.py --ion 6Li --lumi 10 --outdir out
```

Outputs: x–Q² rate maps and per-bin statistical-precision maps for the
three flagship observables, plus `summary_<ion>.txt`.

## Modules

| module | content |
|---|---|
| `beams.py` | species table (d, ³He, ⁶Li, ⁷Li), rigidity-scaled top energies, energy-scan configs |
| `kinematics.py` | DIS variables, scattered-electron lab kinematics, acceptance masks, log grids |
| `structure.py` | **TOY** F2p/F2n (HERA-anchored ±30%), nuclear F2A builder with EMC-ratio hook, NC cross section |
| `polarized.py` | **TOY** g1p/g1n, polarized-EMC ratio hook, toy b1 (HERMES-like), toy Δ(x,Q²) scenarios |
| `asymmetries.py` | spin-1 master formula pieces: A∥ (g1), A_zz (b1), A_cos2φ (Δ); statistical error estimators |
| `fom.py` | luminosity scenarios → events/bin → δ(observable) maps + significance |

## Big caveats (by design — see plans/02, steps 1.3–1.5)

1. Every physics input marked **TOY** must be replaced before quoting
   results: LHAPDF/JAM/DSSV grids for F2/g1, digitized Cloët–Bentz–Thomas
   curves for the polarized-EMC ratio, theory predictions for b1(⁶Li) and
   Δ. The interfaces are built so each swap is one line.
2. Ion in-ring polarizations (`Scenario.pol_ion_*`) are placeholders;
   source targets are P_z ≥ 0.90, P_zz ≥ 0.80 (ECRP proposal Table 1) and
   ring transport survival is an open question (plans/04_open_questions.md).
3. Li top energies are rigidity estimates (275·Z/A GeV/u); confirm with
   BNL C-AD / EPIOS.
4. Statistical errors only; no detector smearing beyond η/E′ cuts of the
   scattered electron (Phase-1 step 1.6 adds eic-smear/Delphes; spectator
   tagging acceptance comes from BeAGLE + far-forward maps in step 1.7).
