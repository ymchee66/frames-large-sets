import sys, time
from math import comb
from ortools.sat.python import cp_model
from km import Grp
from core import lam_star
def partition(A, N, lam, sizes=None, tl=120, workers=2, seed=0):
    m = cp_model.CpModel(); n = len(A[0])
    x = [[m.NewBoolVar(f"x{j}_{c}") for c in range(N)] for j in range(n)]
    for j in range(n): m.AddExactlyOne(x[j])
    for r in A:
        for c in range(N): m.Add(sum(r[j] * x[j][c] for j in range(n) if r[j]) == lam)
    m.Add(x[0][0] == 1)
    s = cp_model.CpSolver(); s.parameters.max_time_in_seconds = tl; s.parameters.num_workers = workers; s.parameters.random_seed = seed
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return "SAT", [[j for j in range(n) if s.Value(x[j][c])] for c in range(N)]
    return ("UNSAT" if st == cp_model.INFEASIBLE else "UNKNOWN"), None
def run(name, p, t, k, tl=120):
    G = Grp(name, p); v = G.v; lam = lam_star(t, k, v); N = comb(v - t, k - t) // lam
    t0 = time.time(); rows, cols, sizes, A = G.km(t, k)
    st, cls = partition(A, N, lam, sizes, tl)
    print(f"{name}(.,{p}) t={t} k={k} v={v} lam*={lam} N={N} rows={len(rows)} cols={len(cols)} -> {st}  [{time.time()-t0:.1f}s]", flush=True)
    return st, cls, rows, cols, sizes, A
if __name__ == "__main__":
    name, p, t, k = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    run(name, p, t, k, float(sys.argv[5]) if len(sys.argv) > 5 else 120)
