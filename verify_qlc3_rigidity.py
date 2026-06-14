# -*- coding: utf-8 -*-
"""Sanity check of the QLC(3) rigidity theorem: alternating projections on the
QLC(3) variety should only ever converge to classical solutions (cardinality 3)."""
import numpy as np

rng = np.random.default_rng(3)

def onb_polar(M):
    U, _, Vh = np.linalg.svd(M)
    return U @ Vh

def cardinality(V, tol=1e-6):
    N = len(V)
    G = np.abs(V @ V.conj().T)
    parent = list(range(N))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for i in range(N):
        for j in range(i + 1, N):
            if G[i, j] > 1 - tol:
                ri, rj = find(i), find(j)
                if ri != rj:
                    parent[ri] = rj
    return len({find(i) for i in range(N)})

v = 3
cards = []
nconv = 0
for run in range(60):
    X = rng.normal(size=(v, v, v, v)) + 1j * rng.normal(size=(v, v, v, v))
    X /= np.linalg.norm(X, axis=3, keepdims=True)
    for it in range(2500):
        for a in range(v):
            for b in range(v):
                X[:, a, b] = onb_polar(X[:, a, b].T).T
                X[a, :, b] = onb_polar(X[a, :, b].T).T
                X[a, b, :] = onb_polar(X[a, b, :].T).T
    res = 0.0
    for a in range(v):
        for b in range(v):
            for M in (X[:, a, b].T, X[a, :, b].T, X[a, b, :].T):
                res = max(res, np.abs(M.conj().T @ M - np.eye(v)).max())
    if res < 1e-9:
        nconv += 1
        cards.append(cardinality(X.reshape(v ** 3, v)))

print(f"converged: {nconv}/60, cardinalities observed: {sorted(set(cards))}")
