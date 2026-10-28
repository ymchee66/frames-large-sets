"""Fast necessary conditions for Theorem 4.1 data (H; K_1..K_m with indices).
N1: every H-orbit of k-sets contains an orbit of some K_i with fusion weight 1.
N4: the H-orbits split into classes U_i, U_i an H-invariant t-(v,k,[H:K_i]*lam) design using only N1-allowed orbits."""
from itertools import combinations
from ortools.sat.python import cp_model
from fix1 import orbits
def n1_n4(v,t,k,lam,Hgens,Ks,idxs,tl=20):
    kH,kHr=orbits(v,k,Hgens); Hsize={}
    for S,d in kH.items(): Hsize[d]=Hsize.get(d,0)+1
    allowed=[set() for _ in Ks]
    for i,K in enumerate(Ks):
        kK,kKr=orbits(v,k,K); Ks_={}
        for S,j in kK.items(): Ks_[j]=Ks_.get(j,0)+1
        for j,S in enumerate(kKr):
            d=kH[S]
            if idxs[i]*Ks_[j]==Hsize[d]: allowed[i].add(d)
    for d in range(len(kHr)):
        if not any(d in a for a in allowed): return "N1"
    tH,tHr=orbits(v,t,Hgens)
    A=[]
    for T in tHr:
        row={}
        for R in combinations([z for z in range(v) if z not in T],k-t):
            d=kH[tuple(sorted(T+R))]; row[d]=row.get(d,0)+1
        A.append(row)
    m=cp_model.CpModel(); x={}
    for i in range(len(Ks)):
        for d in allowed[i]: x[i,d]=m.NewBoolVar("")
    for d in range(len(kHr)): m.AddExactlyOne([x[i,d] for i in range(len(Ks)) if (i,d) in x])
    for i in range(len(Ks)):
        for row in A: m.Add(sum(a*x[i,d] for d,a in row.items() if (i,d) in x)==idxs[i]*lam)
    # symmetry: identical (K,idx) types interchangeable -> order them
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=1
    st=s.Solve(m)
    return {cp_model.INFEASIBLE:"N4",cp_model.OPTIMAL:"pass",cp_model.FEASIBLE:"pass"}.get(st,"N4?")
