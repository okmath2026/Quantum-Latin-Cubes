# -*- coding: utf-8 -*-
"""
Compact deterministic certificate for the QLS(5), c=8 classification.

The certificate is intentionally small: aggregate counts plus a SHA-256 digest
over every canonical pattern and its induced support/orthogonality obstruction.
"""

from __future__ import annotations

from collections import Counter
import hashlib
import json
import os

from c8_pattern_search import (
    N,
    OLD,
    NEW,
    canonicalize_new_labels,
    obstruction_type,
    update_constraints,
)


OUTDIR = os.path.dirname(os.path.abspath(__file__))


def support_tuple(zero_sets):
    return tuple(tuple(i for i in OLD if i not in zero_sets[k]) for k in NEW)


def main():
    grid = [[None] * N for _ in range(N)]
    grid[0] = list(OLD)
    row_used = [set(grid[0])] + [set() for _ in range(1, N)]
    col_used = [{c} for c in OLD]
    zero_sets = {k: set() for k in NEW}
    new_edges = set()
    cells = [(r, c) for r in range(1, N) for c in range(N)]

    seen = set()
    records = []
    counts = Counter()
    support_dims = Counter()
    edge_counts = Counter()
    failures = []

    def record(zero_sets, new_edges):
        key = canonicalize_new_labels(grid)
        if key in seen:
            return
        seen.add(key)
        supports = support_tuple(zero_sets)
        edges = tuple(sorted(new_edges))
        obstruction = obstruction_type(zero_sets, new_edges)
        if obstruction is None:
            failures.append((key, supports, edges))
        counts[obstruction] += 1
        support_dims[tuple(sorted(len(s) for s in supports))] += 1
        edge_counts[len(edges)] += 1
        records.append({
            "pattern": key,
            "supports": supports,
            "edges": edges,
            "obstruction": obstruction,
        })

    def backtrack(idx, zero_sets, new_edges):
        if idx == len(cells):
            if all(any(k in row_used[r] for r in range(1, N)) for k in NEW):
                record(zero_sets, new_edges)
            return

        r, c = cells[idx]
        for val in range(N + 3):
            if val in row_used[r] or val in col_used[c]:
                continue
            if val == c:
                continue
            grid[r][c] = val
            row_used[r].add(val)
            col_used[c].add(val)
            updated = update_constraints(grid, r, c, val, zero_sets, new_edges)
            if updated is not None:
                backtrack(idx + 1, updated[0], updated[1])
            col_used[c].remove(val)
            row_used[r].remove(val)
            grid[r][c] = None

    backtrack(0, zero_sets, new_edges)

    records.sort(key=lambda x: (x["pattern"], x["supports"], x["edges"], str(x["obstruction"])))
    digest = hashlib.sha256()
    for rec in records:
        digest.update(repr(rec).encode("utf-8"))
        digest.update(b"\n")

    cert = {
        "statement": "All normalized QLS(5) cardinality-8 label patterns are obstructed.",
        "normalization": "first row is e0,e1,e2,e3,e4; new projective classes are 5,6,7",
        "unique_patterns": len(records),
        "unobstructed_patterns": len(failures),
        "sha256": digest.hexdigest(),
        "obstruction_counts": {str(k): v for k, v in sorted(counts.items(), key=lambda x: str(x[0]))},
        "support_dimension_counts": {repr(k): v for k, v in sorted(support_dims.items())},
        "new_edge_counts": {str(k): v for k, v in sorted(edge_counts.items())},
    }
    path = os.path.join(OUTDIR, "c8_classification_certificate.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(cert, f, ensure_ascii=False, indent=2)

    print(json.dumps(cert, ensure_ascii=False, indent=2))
    print("wrote", path)
    if failures:
        raise SystemExit("unobstructed patterns found")


if __name__ == "__main__":
    main()
