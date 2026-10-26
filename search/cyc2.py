"""Theorem 3.1 with an arbitrary group G (generators) and shifter sigma on {0..v-1}; brute-force orbits."""
import sys, json, time
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model
from core import *
def cyc_search(v, t, k, N, gens, sigma, tl=300, workers=2, log=True):
    lam = comb(v - t, k - t) // N
    allg = gens + [sigma]
    # orbits of G on k-sets and t-sets
    def orbs(kk):
        idx = {}; reps = []
        for S in combinations(range(v), kk):
            if S in idx: continue
            o = orbit(S, gens) if gens else {S}
            for T in o: idx[T] = len(reps)
            reps.append(S)
        return idx, reps
    kidx, kreps = orbs(k); tidx, treps = orbs(t)
    sig = [kidx[img(S, sigma)] for S in kreps]
    seen = [False] * len(kreps); cycles = []
    for j in range(len(kreps)):
        if not seen[j]:
            c = []; i = j
            while not seen[i]: seen[i] = True; c.append(i); i = sig[i]
            cycles.append(c)
    if any(len(c) != N for c in cycles): return "SHORT", None, len(treps), len(kreps)
    A = [dict() for _ in treps]
    for i, T in enumerate(treps):
        rest = [x for x in range(v) if x not in T]
        for R in combinations(rest, k - t):
            j = kidx[tuple(sorted(T + R))]; A[i][j] = A[i].get(j, 0) + 1
    m = cp_model.CpModel(); x = [m.NewBoolVar(str(j)) for j in range(len(kreps))]
    for c in cycles: m.AddExactlyOne(x[j] for j in c)
    for r in A: m.Add(sum(a * x[j] for j, a in r.items()) == lam)
    m.Add(x[cycles[0][0]] == 1)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "SAT", [list(kreps[j]) for j in range(len(kreps)) if s.Value(x[j])], len(treps), len(kreps)
    return ("UNSAT" if st == cp_model.INFEASIBLE else "UNKNOWN"), None, len(treps), len(kreps)
def shifter(v, N, c):
    """c cycles of length N on points 0..cN-1, the remaining points fixed."""
    return cyc(v, [tuple(range(i * N, (i + 1) * N)) for i in range(c)])
def verify(v, t, k, N, gens, sigma, reps):
    B = set()
    for r in reps: B |= (orbit(tuple(r), gens) if gens else {tuple(sorted(r))})
    lam = comb(v - t, k - t) // N
    ok = is_design(B, v, t) == (True, lam)
    allb = set(B); D = B
    for i in range(1, N):
        D = {img(b, sigma) for b in D}
        if allb & D: ok = False
        allb |= D
    return ok and len(allb) == comb(v, k)
if __name__ == "__main__":
    print("candidates (t,t+1,v): lam*, N")
    for t in (3, 4):
        for v in range(t + 4, 34):
            ls = lam_star(t, t + 1, v); N = (v - t) // ls
            print(t, t + 1, v, ls, N)
