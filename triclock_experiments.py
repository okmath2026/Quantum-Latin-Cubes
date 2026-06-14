# -*- coding: utf-8 -*-
"""
Tri-clock construction for maximal-cardinality quantum Latin cubes of prime order.

C[r,j,k]_s = omega^( r*s + j*pi(s) + k*g(s) ) / sqrt(p)

F1: find (pi, g) with {1, id, pi, g} linearly independent over F_p  (p = 5)
F2: verify the cube for p = 5 and p = 7 (all ONB conditions + cardinality p^3)
G : same permutation condition over Z_v for composite v = 4, 6 (remark material)
"""
import itertools
import numpy as np

def lin_indep_Fp(p, pi, g):
    """Check {1, id, pi, g} linearly independent over F_p (no nontrivial
    (a,b,c,d) with a + b*s + c*pi(s) + d*g(s) = 0 for all s)."""
    M = np.array([[1, s, pi[s], g[s]] for s in range(p)], dtype=np.int64) % p
    # rank over F_p by Gaussian elimination
    A = M.copy() % p
    rank, rows, cols = 0, A.shape[0], A.shape[1]
    r = 0
    for c in range(cols):
        piv = None
        for rr in range(r, rows):
            if A[rr, c] % p != 0:
                piv = rr
                break
        if piv is None:
            continue
        A[[r, piv]] = A[[piv, r]]
        inv = pow(int(A[r, c]), p - 2, p)
        A[r] = (A[r] * inv) % p
        for rr in range(rows):
            if rr != r and A[rr, c] % p != 0:
                A[rr] = (A[rr] - A[rr, c] * A[r]) % p
        r += 1
    return r == 4

def no_relation_Zv(v, pi, g):
    """Over Z_v (possibly composite): no (a,b,c) != (0,0,0) and const t with
    a*s + b*pi(s) + c*g(s) == t (mod v) for all s."""
    s = np.arange(v)
    for a in range(v):
        for b in range(v):
            for c in range(v):
                if a == b == c == 0:
                    continue
                vals = (a * s + b * pi[s] + c * g[s]) % v
                if np.all(vals == vals[0]):
                    return False
    return True

# ----------------------------------------------------------------------
print("=" * 72)
print("F1. search (pi, g) with {1,id,pi,g} linearly independent over F_5")
print("=" * 72)
p = 5
found = []
perms = list(itertools.permutations(range(p)))
pi0 = None
for pi in perms:
    # pi non-affine and pi(0)=0, pi(1)=1 to normalize (not required, just tidy)
    if pi[0] == 0 and pi[1] == 1 and pi != tuple(range(p)):
        for g in perms:
            if lin_indep_Fp(p, np.array(pi), np.array(g)):
                found.append((pi, g))
        if found:
            pi0 = pi
            break
print(f"pi = {pi0},  number of valid g for this pi: {len(found)}")
print(f"first few valid g: {[f[1] for f in found[:5]]}")

# total count over all pairs (sanity of the counting argument)
cnt = 0
for pi in perms:
    for g in perms:
        if lin_indep_Fp(p, np.array(pi), np.array(g)):
            cnt += 1
print(f"total valid (pi,g) pairs over S_5 x S_5: {cnt} / {len(perms)**2}")

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("F2. verify the tri-clock cube for p = 5 and p = 7")
print("=" * 72)

def build_cube(p, pi, g):
    w = np.exp(2j * np.pi / p)
    s = np.arange(p)
    C = np.zeros((p, p, p, p), complex)
    for r in range(p):
        for j in range(p):
            for k in range(p):
                C[r, j, k] = w ** ((r * s + j * pi[s] + k * g[s]) % p) / np.sqrt(p)
    return C

def verify_cube(C, p):
    res = 0.0
    for j in range(p):
        for k in range(p):
            M = C[:, j, k].T
            res = max(res, np.abs(M.conj().T @ M - np.eye(p)).max())
    for r in range(p):
        for k in range(p):
            M = C[r, :, k].T
            res = max(res, np.abs(M.conj().T @ M - np.eye(p)).max())
    for r in range(p):
        for j in range(p):
            M = C[r, j, :].T
            res = max(res, np.abs(M.conj().T @ M - np.eye(p)).max())
    V = C.reshape(p ** 3, p)
    G = np.abs(V @ V.conj().T)
    ncoll = int(((G > 1 - 1e-7).sum() - p ** 3) // 2)
    return res, p ** 3 - 0 if ncoll == 0 else None, ncoll

pi5, g5 = map(np.array, found[0])
C5 = build_cube(5, pi5, g5)
res, _, ncoll = verify_cube(C5, 5)
print(f"p=5: pi={tuple(pi5)}, g={tuple(g5)}")
print(f"p=5: max ONB residual = {res:.2e}, projective coincidences = {ncoll}"
      f"  -> cardinality {125 - ncoll if ncoll == 0 else 'NOT maximal'}")

# p = 7: pi non-affine, g found by counting argument (search quickly)
p = 7
rng = np.random.default_rng(7)
pi7 = np.array([0, 1, 2, 3, 4, 6, 5])   # non-affine (affine maps with pi(0)=0,pi(1)=1 is id)
g7 = None
for trial in range(100000):
    cand = rng.permutation(7)
    if lin_indep_Fp(7, pi7, cand):
        g7 = cand
        break
print(f"p=7: pi={tuple(pi7)}, g={tuple(g7)}")
C7 = build_cube(7, pi7, g7)
res, _, ncoll = verify_cube(C7, 7)
print(f"p=7: max ONB residual = {res:.2e}, projective coincidences = {ncoll}"
      f"  -> cardinality {343 - ncoll if ncoll == 0 else 'NOT maximal'}")

# ----------------------------------------------------------------------
print()
print("=" * 72)
print("G. composite orders: does the permutation condition hold over Z_v?")
print("=" * 72)
for v in (4, 6):
    perms_v = list(itertools.permutations(range(v)))
    ok = None
    tested = 0
    for pi in perms_v:
        for g in perms_v:
            tested += 1
            if no_relation_Zv(v, np.array(pi), np.array(g)):
                ok = (pi, g)
                break
        if ok:
            break
    print(f"v={v}: {'FOUND ' + str(ok) if ok else 'NO pair (pi,g) exists'}"
          f"   (tested {tested} pairs)")
    if ok:
        Cv = build_cube(v, np.array(ok[0]), np.array(ok[1]))
        res, _, ncoll = verify_cube(Cv, v)
        print(f"      cube check: residual {res:.2e}, coincidences {ncoll}, "
              f"cardinality {v**3 - ncoll if ncoll == 0 else 'NOT maximal'}")
