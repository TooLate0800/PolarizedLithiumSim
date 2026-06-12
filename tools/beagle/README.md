# Running BeAGLE locally — status & instructions

Goal: build `../../../BeAGLE` (eic/BeAGLE, master) on this machine for the
e+⁶Li/⁷Li breakup & tagging study (plans/02 step 1.5).

## Dependency status (built 2026-06-12 into `~/Projects/eic/beagle_deps/install`)

| dependency | status | notes |
|---|---|---|
| gfortran 11.4 | ✅ system | legacy code needs `-std=legacy` (handled by build script) |
| LHAPDF 5.9.1 | ✅ built | includes `eps09.f`; grids installed: cteq6ll.LHpdf (=CTEQ6L1, set 10042), EPS09LOR_{4,6,12}.LHgrid — **A=6 is a valid EPS09 nucleus**, and BeAGLE's `anear.f` maps both ⁶Li and ⁷Li → A=6 |
| RAPGAP 3.302 libs | ✅ built | `librapgap33.a` (sfecfe.o removed), `libar4.a`, `libbases.a` copied into `BeAGLE/RAPGAP-3.302/lib/` |
| PYTHIA 6.4.28 | ✅ built | standalone `libpythia6.a` (for RAPGAP configure); BeAGLE builds its own internally |
| CERNLIB 2024 (free) | ✅ built | mathlib/kernlib/packlib static, graphics/PAW tree patched out (no Motif on this box) |
| **FLUKA** | ❌ **needs your action** | license registration is personal — see below |
| nuclear.bin | ✅ in BeAGLE repo | FLUKA evaporation data file (link into run dir) |

## The one missing piece: FLUKA (~10 minutes of your time)

BeAGLE links `libflukahp.a` and compiles against `$FLUPRO/flukapro/`
includes — both come only with the **INFN FLUKA** distribution:

1. Register (free for research) at **https://www.fluka.org** → Download.
2. Download the 64-bit binary release matching **gfortran 11**
   (releases are tagged by gfortran major version, e.g.
   `fluka2024.x-linux-gfor64bit-11_amd64.tar.gz` — pick the gfortran-11
   variant; BeAGLE was developed against the older `fluka2011.2x` line,
   so if the final link complains we install that instead).
3. Unpack and set:
   ```bash
   mkdir -p ~/Projects/eic/beagle_deps/fluka && cd ~/Projects/eic/beagle_deps/fluka
   tar xzf ~/Downloads/fluka20XX...tar.gz
   export FLUPRO=$PWD            # dir containing libflukahp.a and flukapro/
   ```
4. Then: `cd ~/Projects/eic/pol_li/tools/beagle && ./build_beagle.sh`

## Build & run

```bash
./build_beagle.sh            # needs $FLUPRO; builds ./BeAGLE in the BeAGLE repo
source env.sh                # BEAGLESYS, LHAPATH, LD_LIBRARY_PATH
mkdir -p ~/scratch/beagle_test && cd ~/scratch/beagle_test
ln -sf $BEAGLESYS/nuclear.bin .
cp $BEAGLESYS/Examples/eAt1dfJn .                       # PYTHIA card
beagle_run() { $BEAGLESYS/BeAGLE < "$1" > "$1.log"; }
beagle_run ~/Projects/eic/pol_li/tools/beagle/cards/eD_18x135_test.inp   # control first!
beagle_run ~/Projects/eic/pol_li/tools/beagle/cards/eLi7_10x117_draft.inp
```

Validation order (plans/02 step 1.5.3): e+d control vs Tu et al. → fragment
yields for A=6,7 → cluster-IA comparison. Output: PYTHIA-style text on
fort.92? (per `OUTPUT` card) → eic-smear `BuildTree`/`TreeToHepMC` (run
inside eic-shell) → HepMC3.

## No-FLUKA alternatives (available right now)

- **Official samples via xrootd** (verified reachable from this machine
  through the eic-shell container): e+d (`eH2`) and e+³He EVGEN under
  `root://dtn-eic.jlab.org//volatile/eic/EPIC/EVGEN/DIS/BeAGLE1.03.02-*/`
  — enough to develop the full analysis/conversion chain and do the e+d
  control study before any local BeAGLE run.
- **BNL SDCC / JLab ifarm**: BeAGLE prebuilt; submit the e+Li cards there.

## Caveats for e+Li physics (from plans/02 step 1.5)

A>4 uses the C-12 Fermi-momentum parameterization; Woods–Saxon geometry
(no α+d/α+t clustering) — fragment spectra come from FLUKA statistical
de-excitation. Validate, don't trust; the cluster-IA model in
`fastsim/polli_fastsim/spectator.py` brackets the model uncertainty.
