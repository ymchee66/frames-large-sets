"""Cyclic large sets (Theorem 3.1): G = {x -> a x + b : a in M}, M the multiplier subgroup of order m,
shifter sigma = x -> s x with s of order m*N, acting on GF(p) (v=p) or GF(p)+{inf} (v=p+1, inf fixed)."""
import sys, json, time
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model
from core import lam_star
INF = -1
def primroot(p):
    for g in range(2, p):
        x = 1; n = 0
        while True:
            x = x * g % p; n += 1
            if x == 1: break
        if n == p - 1: return g
class AG:
    def __init__(s, p, m, ext):
        s.p, s.m, s.ext = p, m, ext; g = primroot(p)
        h = pow(g, (p - 1) // m, p); s.M = [pow(h, i, p) for i in range(m)]; s.g = g
        s.v = p + ext; s.pts = ([INF] if ext else []) + list(range(p))
    def canon(s, S):
        p = s.p; fin = [x for x in S if x != INF]; best = None
        for a in s.M:
            for z in fin:
                c = tuple(sorted(a * (x - z) % p for x in fin))
                if best is None or c < best: best = c
        return ((INF,) if len(fin) < len(S) else ()) + best
    def reps(s, k):
        out = {}; p = s.p
        for kk, add in ((k, ()),) + (((k - 1, (INF,)),) if s.ext else ()):
            for R in combinations(range(1, p), kk - 1):
                S = add + (0,) + R; c = s.canon(S)
                out.setdefault(c, 0); out[c] += 1
        return sorted(out)
    def km(s, t, k):
        rows = s.reps(t); cols = s.reps(k); idx = {c: j for j, c in enumerate(cols)}
        A = []
        for T in rows:
            r = [0] * len(cols); rest = [x for x in s.pts if x not in T]
            for R in combinations(rest, k - t): r[idx[s.canon(T + R)]] += 1
            A.append(r)
        return rows, cols, idx, A
def search(p, ext, t, k, N, m, tl=120):
    G = AG(p, m, ext); v = G.v; lam = comb(v - t, k - t) // N
    sg = pow(G.g, (p - 1) // (m * N), p)
    rows, cols, idx, A = G.km(t, k)
    sig = [idx[G.canon(tuple(INF if x == INF else sg * x % p for x in c))] for c in cols]
    seen = [False] * len(cols); cycles = []
    for j in range(len(cols)):
        if not seen[j]:
            cyc = []; i = j
            while not seen[i]: seen[i] = True; cyc.append(i); i = sig[i]
            cycles.append(cyc)
    if any(len(c) != N for c in cycles): return "SHORT", None, len(rows), len(cols)
    mdl = cp_model.CpModel(); x = [mdl.NewBoolVar(str(j)) for j in range(len(cols))]
    for c in cycles: mdl.AddExactlyOne(x[j] for j in c)
    for r in A: mdl.Add(sum(r[j] * x[j] for j in range(len(cols)) if r[j]) == lam)
    mdl.Add(x[cycles[0][0]] == 1)
    sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = tl; sv.parameters.num_workers = 2
    st = sv.Solve(mdl)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "SAT", [list(cols[j]) for j in range(len(cols)) if sv.Value(x[j])], len(rows), len(cols)
    return ("UNSAT" if st == cp_model.INFEASIBLE else "UNKNOWN"), None, len(rows), len(cols)
if __name__ == "__main__":
    from gf import is_prime
    pmax = int(sys.argv[1]); res = []
    for p in [q for q in range(7, pmax + 1) if is_prime(q)]:
        for ext in (0, 1):
            v = p + ext
            for t in (2, 3, 4, 5):
                for k in range(t + 1, v // 2 + 1):
                    ls = lam_star(t, k, v); Nmax = comb(v - t, k - t) // ls
                    for N in [d for d in range(3, p) if (p - 1) % d == 0 and Nmax % d == 0]:
                        for m in sorted([d for d in range(1, p) if ((p - 1) // N) % d == 0], reverse=True):
                            if comb(p - 1, k - 1) * k * m > 3e7 or comb(v, k) / (p * m) > 6000: continue
                            t0 = time.time(); st, sol, nr, nc = search(p, ext, t, k, N, m)
                            lam = comb(v - t, k - t) // N
                            print(f"v={v} t={t} k={k} N={N} lam={lam}{'*' if lam==ls else ''} m={m} rows={nr} cols={nc} {st} [{time.time()-t0:.0f}s]", flush=True)
                            if st == "SAT":
                                res.append(dict(p=p, ext=ext, t=t, k=k, N=N, lam=lam, minimal=(lam == ls), m=m, reps=sol))
                                json.dump(res, open("cyc.json", "w")); break
                            if st == "UNKNOWN": break
