"""Paper's 'first approach', made exact: peel off designs one at a time, finish the last few classes by exact partition."""
import sys, json, time, random
from math import comb
from ortools.sat.python import cp_model
from km import Grp
from core import lam_star
from solve import partition
def one_design(A, cols, lam, seed, sizes=None, share=None):
    m = cp_model.CpModel(); x = {j: m.NewBoolVar(f"x{j}") for j in cols}
    if share:
        for sz, cnt in share.items(): m.Add(sum(x[j] for j in cols if sizes[j] == sz) == cnt)
    for r in A: m.Add(sum(r[j] * x[j] for j in cols if r[j]) == lam)
    rnd = random.Random(seed)
    m.Maximize(sum(rnd.randint(0, 5) * x[j] for j in cols))
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = 20; s.parameters.num_workers = 1; s.parameters.random_seed = seed
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE): return [j for j in cols if s.Value(x[j])]
    return None
def peel(A, N, lam, sizes, tail=6, tl=120, tries=50, seed0=0):
    n = len(A[0])
    from collections import Counter
    cnt = Counter(sizes); share = {s: c // N for s, c in cnt.items()} if all(c % N == 0 for c in cnt.values()) else None
    print('share', share, flush=True)
    for tr in range(tries):
        cols = set(range(n)); classes = []; ok = True
        while N - len(classes) > tail:
            d = one_design(A, sorted(cols), lam, seed0 * 1000 + tr * 100 + len(classes), sizes, share)
            if d is None: ok = False; break
            classes.append(d); cols -= set(d)
        if not ok: continue
        rem = sorted(cols); B = [[r[j] for j in rem] for r in A]
        st, cls = partition(B, N - len(classes), lam, tl=tl, seed=tr)
        if st == "SAT":
            return classes + [[rem[j] for j in c] for c in cls]
    return None
if __name__ == "__main__":
    name, t, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); ps = [int(a) for a in sys.argv[4:]]
    for p in ps:
        G = Grp(name, p); v = G.v; lam = lam_star(t, k, v); N = comb(v - t, k - t) // lam
        t0 = time.time(); rows, cols, sizes, A = G.km(t, k)
        cls = peel(A, N, lam, sizes)
        print(f"{name} p={p} t={t} k={k} lam*={lam} N={N} rows={len(rows)} cols={len(cols)} -> {'SAT' if cls else 'FAILED'} [{time.time()-t0:.0f}s]", flush=True)
        if cls:
            json.dump([[list(cols[j]) for j in c] for c in cls], open(f"ls_{name}_{t}_{k}_{p}.json", "w"))
