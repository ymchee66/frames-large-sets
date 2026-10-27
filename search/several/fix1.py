"""Large sets with one design fixed by H=<G,tau> and the other N-1 designs permuted in one cycle by tau.
G fixes every design; tau normalises G and tau^(N-1) is in G."""
import sys, json, time
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model

def compose(a, b):  # apply a then b (right action)
    return tuple(b[a[i]] for i in range(len(a)))
def img(S, g): return tuple(sorted(g[x] for x in S))
def gen_group(gens, v):
    e = tuple(range(v)); G = {e}; fr = [e]
    while fr:
        nw = []
        for g in fr:
            for h in gens:
                gh = compose(g, h)
                if gh not in G: G.add(gh); nw.append(gh)
        fr = nw
    return G
def orbits(v, kk, gens):
    idx = {}; reps = []
    for S in combinations(range(v), kk):
        if S in idx: continue
        o = {S}; fr = [S]
        while fr:
            nw = []
            for T in fr:
                for g in gens:
                    U = img(T, g)
                    if U not in o: o.add(U); nw.append(U)
            fr = nw
        for T in o: idx[T] = len(reps)
        reps.append(S)
    return idx, reps

def search(v, t, k, N, gens, tau, tl=120, workers=4, seed=0):
    lam = comb(v - t, k - t) // N
    Hg = gens + [tau]
    kG, kGr = orbits(v, k, gens); tG, tGr = orbits(v, t, gens)
    kH, kHr = orbits(v, k, Hg); tH, tHr = orbits(v, t, Hg)
    # G-orbits inside each H-orbit
    inside = [[] for _ in kHr]
    for j, S in enumerate(kGr): inside[kH[S]].append(j)
    # KM matrices
    def km(treps, kidx):
        A = []
        for T in treps:
            row = {}; rest = [x for x in range(v) if x not in T]
            for R in combinations(rest, k - t):
                j = kidx[tuple(sorted(T + R))]; row[j] = row.get(j, 0) + 1
            A.append(row)
        return A
    AG = km(tGr, kG); AH = km(tHr, kH)
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{j}") for j in range(len(kGr))]   # B1 (G-orbits)
    y = [m.NewBoolVar(f"y{d}") for d in range(len(kHr))]   # B0 (H-orbits)
    for d, ins in enumerate(inside):
        if len(ins) < N - 1:
            m.Add(y[d] == 1)
            for j in ins: m.Add(x[j] == 0)
        else:
            m.Add(y[d] + sum(x[j] for j in ins) == 1)
    for r in AG: m.Add(sum(a * x[j] for j, a in r.items()) == lam)
    for r in AH: m.Add(sum(a * y[d] for d, a in r.items()) == lam)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers; s.parameters.random_seed = seed
    st = s.Solve(m)
    info = dict(Gk=len(kGr), Hk=len(kHr), short=sum(1 for i in inside if len(i) < N - 1))
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        B1 = [kGr[j] for j in range(len(kGr)) if s.Value(x[j])]
        B0 = [kHr[d] for d in range(len(kHr)) if s.Value(y[d])]
        return "SAT", (B0, B1), info
    return ("UNSAT" if st == cp_model.INFEASIBLE else "UNKNOWN"), None, info

def expand(reps, gens, v):
    G = gen_group(gens, v) if gens else {tuple(range(v))}
    return {img(r, g) for r in reps for g in G}

def verify(v, t, k, N, gens, tau, sol):
    B0r, B1r = sol; lam = comb(v - t, k - t) // N
    B0 = expand(B0r, gens + [tau], v); B1 = expand(B1r, gens, v)
    designs = [B0]; D = B1
    for i in range(N - 1):
        designs.append(D); D = {img(b, tau) for b in D}
    if D != B1: return False
    allk = set()
    for D in designs:
        cnt = {}
        for b in D:
            for T in combinations(b, t): cnt[T] = cnt.get(T, 0) + 1
        if len(cnt) != comb(v, t) or any(c != lam for c in cnt.values()): return False
        if allk & D: return False
        allk |= D
    return len(allk) == comb(v, k)

def cycle_perm(v, cycles):
    p = list(range(v))
    for c in cycles:
        for i, a in enumerate(c): p[a] = c[(i + 1) % len(c)]
    return tuple(p)

if __name__ == "__main__":
    # test: STS(9), N=7, tau of type 6^1 2^1 1^1 and 6^1 1^3
    for cyc in ([tuple(range(6)), (6, 7)], [tuple(range(6))]):
        tau = cycle_perm(9, cyc)
        st, sol, info = search(9, 2, 3, 7, [], tau, tl=60)
        print(cyc, st, info, verify(9, 2, 3, 7, [], tau, sol) if sol else None)
