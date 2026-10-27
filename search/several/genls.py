"""Theorem 4.1 solver: H-invariant LS[N](t,k,v) with base designs having stabilisers K_1..K_m."""
import sys, time
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model
from fix1 import compose, img, gen_group, orbits, cycle_perm


def group_order(gens, v):
    return len(gen_group(gens, v)) if gens else 1


def solve(v, t, k, lam, Hgens, Ks, tl=120, workers=2, seed=0, log=False, fix_first=True):
    N = comb(v - t, k - t) // lam
    oH = group_order(Hgens, v)
    idxs = [oH // group_order(K, v) for K in Ks]
    assert sum(idxs) == N, (idxs, N)
    kH, kHr = orbits(v, k, Hgens)
    Hsize = {}
    for S, d in kH.items():
        Hsize[d] = Hsize.get(d, 0) + 1
    m = cp_model.CpModel()
    cover = [[] for _ in kHr]
    data = []
    for i, K in enumerate(Ks):
        kK, kKr = orbits(v, k, K)
        tK, tKr = orbits(v, t, K)
        Ksize = {}
        for S, j in kK.items():
            Ksize[j] = Ksize.get(j, 0) + 1
        x = [m.NewBoolVar(f"x{i}_{j}") for j in range(len(kKr))]
        for j, S in enumerate(kKr):
            d = kH[S]
            w, r = divmod(idxs[i] * Ksize[j], Hsize[d])
            assert r == 0
            if w == 1:
                cover[d].append(x[j])
            else:
                m.Add(x[j] == 0)
        for T in tKr:
            row = {}
            rest = [z for z in range(v) if z not in T]
            for R in combinations(rest, k - t):
                j = kK[tuple(sorted(T + R))]
                row[j] = row.get(j, 0) + 1
            m.Add(sum(a * x[j] for j, a in row.items()) == lam)
        data.append((kKr, x))
    for d in range(len(kHr)):
        if not cover[d]:
            return "INFEASIBLE(cover)", None, dict(H=len(kHr))
        m.AddExactlyOne(cover[d])
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = tl
    s.parameters.num_workers = workers
    s.parameters.random_seed = seed
    s.parameters.log_search_progress = log
    st = s.Solve(m)
    info = dict(Horb=len(kHr), Korb=[len(d[0]) for d in data], idx=idxs, oH=oH)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "SAT", [[kKr[j] for j in range(len(kKr)) if s.Value(x[j])] for kKr, x in data], info
    return ("UNSAT" if st == cp_model.INFEASIBLE else "UNKNOWN"), None, info


def verify(v, t, k, lam, Hgens, Ks, sol):
    Hel = gen_group(Hgens, v)
    designs = []
    for K, reps in zip(Ks, sol):
        Kel = gen_group(K, v) if K else {tuple(range(v))}
        B = frozenset(img(r, g) for r in reps for g in Kel)
        seen = set()
        for h in Hel:
            D = frozenset(img(b, h) for b in B)
            if D not in seen:
                seen.add(D)
                designs.append(D)
    N = comb(v - t, k - t) // lam
    if len(designs) != N:
        return False
    allk = set()
    for D in designs:
        cnt = {}
        for b in D:
            for T in combinations(b, t):
                cnt[T] = cnt.get(T, 0) + 1
        if len(cnt) != comb(v, t) or any(c != lam for c in cnt.values()):
            return False
        if allk & D:
            return False
        allk |= D
    return len(allk) == comb(v, k)
