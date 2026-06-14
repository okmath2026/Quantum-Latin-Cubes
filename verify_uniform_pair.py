# -*- coding: utf-8 -*-
"""Verify the uniform explicit pair: pi = transposition (p-2, p-1),
g = 3-cycle (0 1 2), for the tri-clock QLC(p) of maximal cardinality."""
import numpy as np

def rank_mod_p(M, p):
    A = M.copy() % p
    rows, cols = A.shape
    r = 0
    for c in range(cols):
        piv = next((rr for rr in range(r, rows) if A[rr, c] % p), None)
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        A[r] = (A[r] * pow(int(A[r, c]), p - 2, p)) % p
        for rr in range(rows):
            if rr != r and A[rr, c] % p:
                A[rr] = (A[rr] - A[rr, c] * A[r]) % p
        r += 1
    return r

for p in (5, 7, 11, 13, 17):
    pi = np.arange(p); pi[p-2], pi[p-1] = p-1, p-2
    g = np.arange(p); g[0], g[1], g[2] = 1, 2, 0
    M = np.array([[1, s, pi[s], g[s]] for s in range(p)], dtype=np.int64)
    li = rank_mod_p(M, p) == 4

    w = np.exp(2j*np.pi/p); s = np.arange(p)
    C = np.zeros((p, p, p, p), complex)
    for r_ in range(p):
        for j in range(p):
            for k in range(p):
                C[r_, j, k] = w ** ((r_*s + j*pi[s] + k*g[s]) % p) / np.sqrt(p)
    res = 0.0
    for a in range(p):
        for b in range(p):
            for M_ in (C[:, a, b].T, C[a, :, b].T, C[a, b, :].T):
                res = max(res, np.abs(M_.conj().T @ M_ - np.eye(p)).max())
    V = C.reshape(p**3, p)
    # count projective coincidences blockwise to save memory
    ncoll = 0
    for i0 in range(0, p**3, 512):
        blk = V[i0:i0+512]
        G = np.abs(blk @ V.conj().T)
        ncoll += int((G > 1 - 1e-7).sum()) - blk.shape[0]
    ncoll //= 2
    print(f"p={p:2d}: lin.indep={li}, ONB residual={res:.1e}, "
          f"coincidences={ncoll}, cardinality={'MAXIMAL ' + str(p**3) if ncoll == 0 else 'NOT MAXIMAL'}")
