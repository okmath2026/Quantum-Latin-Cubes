# The Maximal Cardinality of Quantum Latin Cubes

This repository contains the manuscript, verification scripts, and selected
machine-readable certificates for the paper

**The Maximal Cardinality of Quantum Latin Cubes: A Complete Solution**.

The main result proves that a quantum Latin cube (QLC) of order `v` with
maximal cardinality `v^3` exists if and only if `v = 1` or `v >= 4`.

## Contents

- `quantum_latin_cubes.tex` - LaTeX source of the manuscript.
- `quantum_latin_cubes.pdf` - compiled manuscript.
- `verify_uniform_pair.py` - checks the prime tri-clock construction for
  selected primes.
- `verify_composite.py` - checks the composite twisted-product construction
  for selected orders.
- `verify_qlc3_rigidity.py` - checks the order-3 rigidity calculation.
- `triclock_experiments.py` - auxiliary tri-clock experiments.
- `qlc_experiments.py` - exploratory QLC spectrum experiments.
- `qls_spectrum/` - selected scripts and certificates for the QLS(5)
  low-cardinality analysis.

## Basic Verification

Run from the repository root:

```bash
python verify_uniform_pair.py
python verify_composite.py
python verify_qlc3_rigidity.py
python qls_spectrum/exact_qls5_c12_construction.py
```

See `REPRODUCIBILITY.md` for the full list of scripts and certificate files
referenced by the manuscript's computational note.

The manuscript can be compiled with a standard LaTeX engine supporting
`amsart`, `amsmath`, `amssymb`, `amsthm`, `geometry`, `booktabs`, and
`url`, and `hyperref`.

## Notes

The QLS(5) low-cardinality nonexistence component is computer-assisted.
The included JSON files are selected certificates referenced by the
manuscript's computational note.
