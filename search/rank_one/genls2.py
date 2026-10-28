"""Theorem 4.1 solver with explicit indices [H:K_i] (no group enumeration) and design-level verification."""
from itertools import combinations
from math import comb
from ortools.sat.python import cp_model
from fix1 import img, orbits
def solve(v,t,k,lam,Hgens,Ks,idxs,tl=120,workers=2,seed=0):
    kH,kHr=orbits(v,k,Hgens); Hsize={}
    for S,d in kH.items(): Hsize[d]=Hsize.get(d,0)+1
    m=cp_model.CpModel(); cover=[[] for _ in kHr]; data=[]
    for i,K in enumerate(Ks):
        kK,kKr=orbits(v,k,K); tK,tKr=orbits(v,t,K); Ksize={}
        for S,j in kK.items(): Ksize[j]=Ksize.get(j,0)+1
        x=[m.NewBoolVar("") for _ in kKr]
        for j,S in enumerate(kKr):
            d=kH[S]; w,r=divmod(idxs[i]*Ksize[j],Hsize[d])
            if r: raise ValueError("index/orbit mismatch: wrong index for K_%d"%i)
            if w==1: cover[d].append(x[j])
            else: m.Add(x[j]==0)
        for T in tKr:
            row={}
            for R in combinations([z for z in range(v) if z not in T],k-t):
                j=kK[tuple(sorted(T+R))]; row[j]=row.get(j,0)+1
            m.Add(sum(a*x[j] for j,a in row.items())==lam)
        data.append((kKr,x))
    for d in range(len(kHr)):
        if not cover[d]: return "INFEASIBLE(cover)",None,dict(Horb=len(kHr))
        m.AddExactlyOne(cover[d])
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=workers; s.parameters.random_seed=seed
    st=s.Solve(m); info=dict(Horb=len(kHr),Korb=[len(d[0]) for d in data])
    if st in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        return "SAT",[[kKr[j] for j in range(len(kKr)) if s.Value(x[j])] for kKr,x in data],info
    return ("UNSAT" if st==cp_model.INFEASIBLE else "UNKNOWN"),None,info
def expand(v,gens,reps):
    """orbit of the base blocks under <gens>"""
    B=set(reps); fr=list(reps)
    while fr:
        nw=[]
        for b in fr:
            for g in gens:
                c=img(b,g)
                if c not in B: B.add(c); nw.append(c)
        fr=nw
    return frozenset(B)
def verify(v,t,k,lam,Hgens,Ks,sol,N):
    designs=set()
    for K,reps in zip(Ks,sol):
        D0=expand(v,K,reps); orb={D0}; fr=[D0]
        while fr:
            nw=[]
            for D in fr:
                for g in Hgens:
                    E=frozenset(img(b,g) for b in D)
                    if E not in orb: orb.add(E); nw.append(E)
            fr=nw
        designs|=orb
    if len(designs)!=N: return False
    allk=set()
    for D in designs:
        c={}
        for b in D:
            for T in combinations(b,t): c[T]=c.get(T,0)+1
        if len(c)!=comb(v,t) or set(c.values())!={lam} or allk&D: return False
        allk|=D
    return len(allk)==comb(v,k)
