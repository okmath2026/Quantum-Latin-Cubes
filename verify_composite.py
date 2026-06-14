# -*- coding: utf-8 -*-
"""Verify Theorem (composite case): twisted-layer cube over the flat product
square, with explicit sqrt(2)-phases and the unit-double-difference assignment g.
Tested for v = 6 (2x3), 9 (3x3), 10 (2x5)."""
import numpy as np

def Fmat(n):
    w = np.exp(2j * np.pi / n)
    s = np.arange(n)
    return w ** np.outer(s, s) / np.sqrt(n)

def flat_basis(n, i):
    """B(phi^(i)) with phi^(i) = (sqrt2 * i, 0, ..., 0)."""
    D = np.ones(n, complex)
    D[0] = np.exp(1j * np.sqrt(2) * i)
    return D[:, None] * Fmat(n)        # columns are the basis vectors

def build(m, n):
    v = m * n
    U = np.zeros((m, n, n), complex)   # u_{i,a} in C^n
    for i in range(m):
        B = flat_basis(n, i)
        for a in range(n):
            U[i, a] = B[:, a]
    V = np.zeros((n, m, m), complex)   # v_{j,b} in C^m
    for j in range(n):
        B = flat_basis(m, j)
        for b in range(m):
            V[j, b] = B[:, b]
    # product square W (ZZLC Lemma 2.2): block (i,j), inner (k,l)
    W = np.zeros((v, v, v), complex)
    for i in range(m):
        for k in range(n):
            for j in range(n):
                for l in range(m):
                    W[i * n + k, j * m + l] = np.kron(U[i, (j + k) % n], V[j, (i + l) % m])
    # assignment g on Z_n x Z_m with g(0,0)=0, g(1,0)=1, g(0,1)=2, g(1,1)=4
    cells = [(0, 0), (1, 0), (0, 1), (1, 1)]
    vals = [0, 1, 2, 4]
    rest_cells = [(s, t) for s in range(n) for t in range(m) if (s, t) not in cells]
    rest_vals = [x for x in range(v) if x not in vals]
    g = {}
    for c, val in zip(cells + rest_cells, vals + rest_vals):
        g[c] = val
    lam = np.zeros(v, complex)
    for (s, t), val in g.items():
        lam[s * m + t] = np.exp(2j * np.pi * val / v)   # kron index = s*m + t
    # cube
    C = np.zeros((v, v, v, v), complex)
    for kappa in range(v):
        C[:, :, kappa] = W * (lam ** kappa)[None, None, :]
    return W, C, v

for (m, n) in ((2, 3), (3, 3), (2, 5)):
    W, C, v = build(m, n)
    res = 0.0
    for r in range(v):
        M = W[r].T; res = max(res, np.abs(M.conj().T @ M - np.eye(v)).max())
        M = W[:, r].T; res = max(res, np.abs(M.conj().T @ M - np.eye(v)).max())
    Wf = W.reshape(v * v, v)
    G = np.abs(Wf @ Wf.conj().T)
    cW = v * v - (int((G > 1 - 1e-7).sum()) - v * v) // 2
    resC = 0.0
    for a in range(v):
        for b in range(v):
            for M in (C[:, a, b].T, C[a, :, b].T, C[a, b, :].T):
                resC = max(resC, np.abs(M.conj().T @ M - np.eye(v)).max())
    Cf = C.reshape(v ** 3, v)
    ncoll = 0
    for i0 in range(0, v ** 3, 512):
        blk = Cf[i0:i0 + 512]
        G = np.abs(blk @ Cf.conj().T)
        ncoll += int((G > 1 - 1e-7).sum()) - blk.shape[0]
    ncoll //= 2
    print(f"v={v} (m={m},n={n}): W residual {res:.1e}, card(W)={cW} (target {v*v}); "
          f"cube residual {resC:.1e}, coincidences {ncoll}, "
          f"cardinality {'MAXIMAL ' + str(v**3) if ncoll == 0 else 'NOT MAXIMAL'}")
