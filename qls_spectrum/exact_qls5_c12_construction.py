# -*- coding: utf-8 -*-
"""
Exact integer-vector construction of a QLS(5) with cardinality 12.

Vectors are stored unnormalized.  Orthogonality and projective cardinality are
checked exactly using integer dot products.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
import json
import os


OUTDIR = os.path.dirname(os.path.abspath(__file__))
N = 5


def e(i):
    return tuple(1 if j == i else 0 for j in range(N))


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def projective_key(v):
    first = next(x for x in v if x)
    return tuple(Fraction(x, first) for x in v)


def verify(W):
    for r in range(N):
        for c1, c2 in combinations(range(N), 2):
            if dot(W[r][c1], W[r][c2]) != 0:
                raise AssertionError(("row", r, c1, c2, dot(W[r][c1], W[r][c2])))
    for c in range(N):
        for r1, r2 in combinations(range(N), 2):
            if dot(W[r1][c], W[r2][c]) != 0:
                raise AssertionError(("col", c, r1, r2, dot(W[r1][c], W[r2][c])))


def cardinality(W):
    return len({projective_key(v) for row in W for v in row})


def main():
    # Work inside span(e1,e3,e4).  In coordinates (u,v,w)=(e1,e3,e4):
    # c=(1,0,1), evec=(2,1,0), a=(-1,2,1);
    # b=(1,1,-1), d=(-1,2,-5), f=(1,-2,0), g=(1,0,-1).
    a = (0, -1, 0, 2, 1)
    b = (0, 1, 0, 1, -1)
    c = (0, 1, 0, 0, 1)
    d = (0, -1, 0, 2, -5)
    h = (0, 2, 0, 1, 0)   # called e in the derivation; avoid shadowing basis e()
    f = (0, 1, 0, -2, 0)
    g = (0, 1, 0, 0, -1)

    W = [
        [e(0), e(1), e(2), e(3), e(4)],
        [a,    e(0), b,    c,    e(2)],
        [d,    e(2), a,    e(0), h],
        [h,    e(4), e(0), e(2), f],
        [e(2), e(3), c,    g,    e(0)],
    ]

    verify(W)
    card = cardinality(W)
    if card != 12:
        raise AssertionError(card)

    out = {
        "cardinality": card,
        "matrix": [[list(v) for v in row] for row in W],
        "labels": [
            ["e0", "e1", "e2", "e3", "e4"],
            ["a", "e0", "b", "c", "e2"],
            ["d", "e2", "a", "e0", "h"],
            ["h", "e4", "e0", "e2", "f"],
            ["e2", "e3", "c", "g", "e0"],
        ],
        "vectors": {
            "a": list(a),
            "b": list(b),
            "c": list(c),
            "d": list(d),
            "h": list(h),
            "f": list(f),
            "g": list(g),
        },
    }
    path = os.path.join(OUTDIR, "exact_qls5_c12_construction.json")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=2)

    print("verified exact QLS(5) cardinality", card)
    print("wrote", path)


if __name__ == "__main__":
    main()
