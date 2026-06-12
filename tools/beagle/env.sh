# Runtime environment for locally built BeAGLE. Usage: source env.sh
export BEAGLESYS="${BEAGLESYS:-$HOME/Projects/eic/BeAGLE}"
export DEPS="${DEPS:-$HOME/Projects/eic/beagle_deps/install}"
export LHAPATH="$DEPS/share/lhapdf/PDFsets"
export LD_LIBRARY_PATH="$DEPS/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
# FLUPRO must point at the INFN FLUKA install (user-licensed); needed at
# build time and for locating auxiliary data at run time.
[ -n "${FLUPRO:-}" ] || echo "note: FLUPRO not set (only needed to (re)build)"
echo "BEAGLESYS=$BEAGLESYS"
echo "LHAPATH=$LHAPATH"
