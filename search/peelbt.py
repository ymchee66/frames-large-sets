"""Peel designs one at a time with backtracking; the last class is automatically a design."""
import sys, json, time, random
from math import comb
from collections import Counter
from km import Grp, GrpInf
from core import lam_star
from peel import one_design
def peel_bt(A, N, lam, sizes, back=3, log=True):
    n = len(A[0]); cnt = Counter(sizes)
    share = {s: c // N for s, c in cnt.items()} if all(c % N == 0 for c in cnt.values()) else None
    cols = set(range(n)); classes = []; seed = 0; fails = 0
    while len(classes) < N - 1:
        seed += 1
        d = one_design(A, sorted(cols), lam, seed, sizes, share)
        if d is None:
            fails += 1
            for _ in range(min(len(classes), back + fails // 5)):
                cols |= set(classes.pop())
            if log: print(f"  stuck; back to {len(classes)} classes (fail {fails})", flush=True)
            if fails > 200: return None
            continue
        classes.append(d); cols -= set(d)
    last = sorted(cols)
    assert all(sum(r[j] for j in last) == lam for r in A)
    return classes + [last]
if __name__ == "__main__":
    name, t, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    for p in [int(a) for a in sys.argv[4:]]:
        G = GrpInf(name, p) if name.endswith("i") else Grp(name, p)
        v = G.v; lam = lam_star(t, k, v); N = comb(v - t, k - t) // lam
        t0 = time.time(); rows, cols, sizes, A = G.km(t, k)
        print(f"{name} p={p} {t}-({v},{k},{lam}) N={N} rows={len(rows)} cols={len(cols)} KM [{time.time()-t0:.0f}s]", flush=True)
        cls = peel_bt(A, N, lam, sizes)
        print(f"{name} p={p} {t}-({v},{k},{lam}) N={N} -> {'SAT' if cls else 'FAILED'} [{time.time()-t0:.0f}s]", flush=True)
        if cls: json.dump([[list(cols[j]) for j in c] for c in cls], open(f"ls_{name}_{t}_{k}_{p}.json", "w"))
