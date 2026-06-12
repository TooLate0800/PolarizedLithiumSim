#!/bin/bash
# Build BeAGLE against the locally built deps (see README.md).
# Requires: FLUPRO pointing at an INFN FLUKA install (libflukahp.a + flukapro/).
# All Makefile paths are overridden on the make command line -- no edits to
# the BeAGLE checkout needed.
set -euo pipefail

BEAGLE_DIR="${BEAGLE_DIR:-$HOME/Projects/eic/BeAGLE}"
DEPS="${DEPS:-$HOME/Projects/eic/beagle_deps/install}"

if [ -z "${FLUPRO:-}" ] || [ ! -f "$FLUPRO/libflukahp.a" ]; then
    echo "ERROR: FLUPRO is unset or \$FLUPRO/libflukahp.a not found."
    echo "       Register at https://www.fluka.org, download the gfortran-11"
    echo "       64-bit release, unpack, and export FLUPRO=<that dir>."
    exit 1
fi
for f in "$DEPS/lib/libLHAPDF.so" "$DEPS/lib/libmathlib.a" \
         "$BEAGLE_DIR/RAPGAP-3.302/lib/librapgap33.a"; do
    [ -e "$f" ] || { echo "ERROR: missing dependency $f (see README.md)"; exit 1; }
done

# Original BeAGLE flags + -std=legacy (gfortran>=10 argument-mismatch) + -w
FFLAGS="-c -g -m64 -fno-inline -fno-automatic -fno-align-commons -std=legacy -w"

cd "$BEAGLE_DIR"
make all \
    FLUKA="$FLUPRO" \
    LIB1="-L$DEPS/lib -lmathlib -lkernlib -lpacklib -ldl -lm" \
    LIB3="-L$DEPS/lib -lLHAPDF" \
    CXXFLAGS="$FFLAGS" \
    "$@"

echo
echo "BeAGLE built: $BEAGLE_DIR/BeAGLE"
echo "Next: source $(dirname "$0")/env.sh and run a control card (README.md)."
