"""Grid-backend tests (skipped automatically if parton/grids are absent)."""

import pathlib
import sys

import numpy as np
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

parton = pytest.importorskip("parton")


def _have(setname):
    try:
        from parton import mkPDF
        mkPDF(setname, 0)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _have("CT18NLO"), reason="CT18NLO grid not installed")
def test_parton_f2_sane():
    from polli_fastsim.structure import PartonF2
    f2 = PartonF2()
    # HERA-anchored magnitudes and n/p < 1 at valence x
    assert 0.8 < f2.f2p(1e-3, 10.0) < 1.3
    assert 0.3 < f2.f2p(0.1, 10.0) < 0.55
    assert f2.f2n(0.5, 10.0) < f2.f2p(0.5, 10.0)
    # vectorized call works on arrays
    vals = f2.f2p(np.array([1e-3, 0.1, 0.5]), np.array([10.0, 10.0, 10.0]))
    assert vals.shape == (3,) and np.all(vals > 0)


@pytest.mark.skipif(not _have("NNPDFpol11_100"),
                    reason="NNPDFpol11 grid not installed")
def test_parton_g1_sane():
    from polli_fastsim.polarized import PartonG1
    g1 = PartonG1()
    # proton g1 positive at mid x; A1p grows with x and stays below 1
    assert g1.g1p(0.1, 10.0) > 0
    f1p = g1.base.f2p(0.3, 10.0) / (2 * 0.3 * 1.18)
    a1p = float(g1.g1p(0.3, 10.0)) / f1p
    assert 0.2 < a1p < 0.9
    # neutron g1 small/negative at low-mid x (isospin swap sanity)
    assert g1.g1n(0.05, 10.0) < g1.g1p(0.05, 10.0)
