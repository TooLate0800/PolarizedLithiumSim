"""Backend selector: TOY parameterizations vs LHAPDF grids (via parton)."""


def get_backends(pdf="toy"):
    """Return dict with 'base' (F2 source), 'g1' (g1 model), 'tag'."""
    if pdf == "grid":
        from .polarized import PartonG1
        from .structure import PartonF2
        base = PartonF2()  # CT18NLO
        return {"base": base, "g1": PartonG1(base=base), "tag": "grid"}
    from .polarized import ToyG1
    from .structure import ToyF2
    base = ToyF2()
    return {"base": base, "g1": ToyG1(base=base), "tag": "toy"}
