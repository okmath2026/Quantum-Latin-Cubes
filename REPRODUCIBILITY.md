# Reproducibility Guide

This file lists the scripts and certificate files supporting the computational
checks in the manuscript.  The repository is intended for code and selected
machine-readable certificates, not for storing the manuscript itself.

## Cube Construction Checks

Run these commands from the repository root.

```bash
python verify_uniform_pair.py
python verify_composite.py
python verify_qlc3_rigidity.py
python triclock_experiments.py
```

- `verify_uniform_pair.py` checks the prime tri-clock construction for
  `p = 5, 7, 11, 13, 17`.
- `verify_composite.py` checks the twisted-product construction for
  `v = 6, 9, 10`.
- `verify_qlc3_rigidity.py` checks the order-3 rigidity calculation.
- `triclock_experiments.py` records the exhaustive order-4 search, the
  order-6 tri-clock pair, and the count of valid pairs in `S_5 x S_5`.

## Order-Five Square-Spectrum Certificates

The selected scripts and certificates for the low-cardinality QLS(5)
analysis are stored in `qls_spectrum/`.

Primary entry points:

```bash
python qls_spectrum/c8_classification_certificate.py
python qls_spectrum/qls5_c9_deterministic_certificate.py
python qls_spectrum/qls5_c10_deterministic_probe.py
python qls_spectrum/qls5_incremental_probe.py
python qls_spectrum/qls_support_obstruction_general.py
python qls_spectrum/exact_qls5_c12_construction.py
```

Selected machine-readable outputs:

- `qls_spectrum/c8_classification_certificate.json`
- `qls_spectrum/qls5_c9_deterministic_certificate.json`
- `qls_spectrum/qls5_c10_deterministic_probe_complete.json`
- `qls_spectrum/qls5_c11_incremental_old_and_new_ok_complete_reclassified.json`
- `qls_spectrum/exact_qls5_c12_construction.json`
