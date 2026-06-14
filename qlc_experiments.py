# -*- coding: utf-8 -*-
"""
Numerical experiments for the paper on maximal-cardinality quantum Latin cubes.

A: explicit verification of the twisted-layer QLC(6) construction (composite case)
B: cyclic phase system for a flat QLS(5) of maximal cardinality (prime case, layered route)
C: general flat QLS(5) search (no pattern assumption)
D: direct QLC(5) feasibility by alternating projections (general, non-layered route)
E: cyclic phase system for flat QLS(7)
"""
import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(20260612)

def Fmat(n):
    w = np.exp(2j * np.pi / n)
    s = np.arange(n)
    return w ** np.outer(s, s) / np.sqrt(n)

def onb_polar(M):
    """Nearest unitary (polar factor) to square matrix M (columns = vectors)."""
    U, _, Vh = np.linalg.svd(M)
    return U @ Vh

def cardinality(vecs, tol=1e-6):
    """Number of projective classes among unit vectors (rows of vecs)."""
    Nv = len(vecs)
    G = np.abs(vecs @ vecs.conj().T)
    same = G > 1 - tol
    # union-find
    parent = list(range(Nv))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a in range(Nv):
        for b in range(a + 1, Nv):
            if same[a, b]:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
    return len({find(a) for a in range(Nv)})

# ----------------------------------------------------------------------
print("=" * 72)
print("A. Twisted-layer QLC(6) with maximal cardinality 216")
print("=" * 72)

n, m = 3, 2          # u-side in C^3, v-side in C^2, order v = mn = 6
F3, F2 = Fmat(3), Fmat(2)

# m = 2 flat ONBs of C^3 (rows of the m x n row-quantum Latin rectangle U)
phis = [rng.uniform(0, 2 * np.pi, 3) for _ in range(m)]
# n = 3 flat ONBs of C^2 (rows of the n x m row-quantum Latin rectangle V)
psis = [rng.uniform(0, 2 * np.pi, 2) for _ in range(n)]

U = np.zeros((m, n, 3), complex)     # u_{i,k} in C^3
for i in range(m):
    B = np.diag(np.exp(1j * phis[i])) @ F3
    for k in range(n):
        U[i, k] = B[:, k]
V = np.zeros((n, m, 2), complex)     # v_{j,l} in C^2
for j in range(n):
    B = np.diag(np.exp(1j * psis[j])) @ F2
    for l in range(m):
        V[j, l] = B[:, l]

# ZZLC Lemma 2.2 product square W: row (i,k), column (j,l)
W = np.zeros((6, 6, 6), complex)
for i in range(m):
    for k in range(n):
        for j in range(n):
            for l in range(m):
                W[i * n + k, j * m + l] = np.kron(U[i, (j + k) % n], V[j, (i + l) % m])

# check W is a QLS(6) of cardinality 36
res = 0.0
for r in range(6):
    Mr = W[r].T          # columns = row entries
    res = max(res, np.abs(Mr.conj().T @ Mr - np.eye(6)).max())
for c in range(6):
    Mc = W[:, c].T
    res = max(res, np.abs(Mc.conj().T @ Mc - np.eye(6)).max())
print(f"W row/col orthonormality residual : {res:.2e}")
print(f"W cardinality (target 36)         : {cardinality(W.reshape(36, 6))}")

# eigenvalue assignment g on basis index e = 2 s + t, s in Z3, t in Z2
g = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 4, (2, 0): 3, (2, 1): 5}
lam = np.zeros(6, complex)
for (s, t), val in g.items():
    lam[2 * s + t] = np.exp(2j * np.pi * val / 6)

# cube: C[r, c, kappa] = diag(lam^kappa) W[r, c]
C = np.zeros((6, 6, 6, 6), complex)
for kappa in range(6):
    C[:, :, kappa] = W * (lam ** kappa)[None, None, :]

res = 0.0
for kappa in range(6):
    for r in range(6):
        Mr = C[r, :, kappa].T
        res = max(res, np.abs(Mr.conj().T @ Mr - np.eye(6)).max())
    for c in range(6):
        Mc = C[:, c, kappa].T
        res = max(res, np.abs(Mc.conj().T @ Mc - np.eye(6)).max())
for r in range(6):
    for c in range(6):
        Mf = C[r, c].T
        res = max(res, np.abs(Mf.conj().T @ Mf - np.eye(6)).max())
print(f"cube row/col/file residual        : {res:.2e}")
print(f"cube cardinality (target 216)     : {cardinality(C.reshape(216, 6))}")

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("B. Cyclic phase system for a flat QLS(5) of maximal cardinality")
print("=" * 72)
# rows r = 0..4 : entry (r,j) = D_{phi_r} f_{j+r},  phi_0 = 0
# column orthogonality <=> for all r<r' :
#   G(phi_{r'} - phi_r, r'-r) = (1/5) sum_s exp(i (phi_{r'}-phi_r)_s) w^{s(r'-r)} = 0

w5 = np.exp(2j * np.pi / 5)
pairs5 = [(r, rp) for r in range(5) for rp in range(r + 1, 5)]

def resid_cyclic5(x):
    phi = np.zeros((5, 5))
    phi[1:] = x.reshape(4, 5)
    out = []
    for (r, rp) in pairs5:
        d = phi[rp] - phi[r]
        s = np.arange(5)
        val = np.sum(np.exp(1j * d) * w5 ** (s * (rp - r))) / 5
        out += [val.real, val.imag]
    return np.array(out)

best = None
NTRY = 4000
for t in range(NTRY):
    x0 = rng.uniform(0, 2 * np.pi, 20)
    sol = least_squares(resid_cyclic5, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if best is None or sol.cost < best.cost:
        best = sol
print(f"tries = {NTRY},  best residual 2-norm = {np.sqrt(2*best.cost):.3e}")
if np.sqrt(2 * best.cost) < 1e-10:
    phi = np.zeros((5, 5)); phi[1:] = best.x.reshape(4, 5)
    Q = np.zeros((5, 5, 5), complex)
    F5 = Fmat(5)
    for r in range(5):
        for j in range(5):
            Q[r, j] = np.exp(1j * phi[r]) * F5[:, (j + r) % 5] * np.sqrt(5) / np.sqrt(5)
    print(f"flat QLS(5) FOUND, cardinality = {cardinality(Q.reshape(25,5))}")
else:
    print("no solution found -> evidence the cyclic flat system is infeasible")

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("C. General flat QLS(5): least-squares over all 125 entry phases")
print("=" * 72)
# entries q[r,c]_s = exp(i Theta[r,c,s]) / sqrt(5); rows+columns must be ONBs.

def resid_flat5(x):
    Th = x.reshape(5, 5, 5)
    Q = np.exp(1j * Th) / np.sqrt(5)
    out = []
    for r in range(5):
        M = Q[r]                    # 5 vectors (rows of M)
        G = M.conj() @ M.T
        for a in range(5):
            for b in range(a + 1, 5):
                out += [G[a, b].real, G[a, b].imag]
    for c in range(5):
        M = Q[:, c]
        G = M.conj() @ M.T
        for a in range(5):
            for b in range(a + 1, 5):
                out += [G[a, b].real, G[a, b].imag]
    return np.array(out)

results = []
NTRY = 250
for t in range(NTRY):
    x0 = rng.uniform(0, 2 * np.pi, 125)
    sol = least_squares(resid_flat5, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    rn = np.sqrt(2 * sol.cost)
    if rn < 1e-10:
        Th = sol.x.reshape(5, 5, 5)
        Q = np.exp(1j * Th) / np.sqrt(5)
        results.append(cardinality(Q.reshape(25, 5)))
hist = {}
for cval in results:
    hist[cval] = hist.get(cval, 0) + 1
print(f"converged runs: {len(results)} / {NTRY}; cardinality histogram: {dict(sorted(hist.items()))}")
print("max flat-QLS(5) cardinality observed:", max(results) if results else None)

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("D. Direct QLC(5) feasibility: alternating projections on 3-axis ONB sets")
print("=" * 72)

def ap_qlc(v=5, iters=4000, runs=40):
    outcomes = []
    for run in range(runs):
        X = rng.normal(size=(v, v, v, v)) + 1j * rng.normal(size=(v, v, v, v))
        X /= np.linalg.norm(X, axis=3, keepdims=True)
        for it in range(iters):
            # axis 0 lines: fix (j,k), vary i
            for j in range(v):
                for k in range(v):
                    X[:, j, k] = onb_polar(X[:, j, k].T).T
            # axis 1 lines
            for i in range(v):
                for k in range(v):
                    X[i, :, k] = onb_polar(X[i, :, k].T).T
            # axis 2 lines
            for i in range(v):
                for j in range(v):
                    X[i, j, :] = onb_polar(X[i, j, :].T).T
        # residual
        res = 0.0
        for j in range(v):
            for k in range(v):
                M = X[:, j, k].T
                res = max(res, np.abs(M.conj().T @ M - np.eye(v)).max())
        for i in range(v):
            for k in range(v):
                M = X[i, :, k].T
                res = max(res, np.abs(M.conj().T @ M - np.eye(v)).max())
        card = cardinality(X.reshape(v**3, v)) if res < 1e-9 else None
        outcomes.append((res, card))
        print(f"  run {run:2d}: residual {res:.2e}  cardinality {card}")
    return outcomes

out = ap_qlc(v=5, iters=3000, runs=12)
feas = [c for (r, c) in out if c is not None]
print(f"converged: {len(feas)}/12, cardinalities: {sorted(set(feas)) if feas else 'NONE'}")

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("E. Cyclic phase system for a flat QLS(7)")
print("=" * 72)
w7 = np.exp(2j * np.pi / 7)
pairs7 = [(r, rp) for r in range(7) for rp in range(r + 1, 7)]

def resid_cyclic7(x):
    phi = np.zeros((7, 7))
    phi[1:] = x.reshape(6, 7)
    out = []
    s = np.arange(7)
    for (r, rp) in pairs7:
        d = phi[rp] - phi[r]
        val = np.sum(np.exp(1j * d) * w7 ** (s * (rp - r))) / 7
        out += [val.real, val.imag]
    return np.array(out)

best7 = None
NTRY = 3000
for t in range(NTRY):
    x0 = rng.uniform(0, 2 * np.pi, 42)
    sol = least_squares(resid_cyclic7, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15)
    if best7 is None or sol.cost < best7.cost:
        best7 = sol
print(f"tries = {NTRY},  best residual 2-norm = {np.sqrt(2*best7.cost):.3e}")
if np.sqrt(2 * best7.cost) < 1e-10:
    phi = np.zeros((7, 7)); phi[1:] = best7.x.reshape(6, 7)
    F7 = Fmat(7)
    Q = np.zeros((7, 7, 7), complex)
    for r in range(7):
        for j in range(7):
            Q[r, j] = np.exp(1j * phi[r]) * F7[:, (j + r) % 7]
    print(f"flat QLS(7) FOUND, cardinality = {cardinality(Q.reshape(49,7))}")
else:
    print("no solution found for p = 7 cyclic system")
