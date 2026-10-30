"""Certify failure of Prop 4.1 for PSL(3,q) frames: find a k-subset S that is a union of orbits of
the monomial group U (triangle stabiliser), which fixes no point and no line of PG(2,q) for q>2."""
import json,sys
from psl3 import groups
def closure(gens,n):
    seen={tuple(range(n))}; fr=list(seen)
    while fr:
        g=fr.pop(); 
        for h in gens:
            x=tuple(h[g[i]] for i in range(n))
            if x not in seen: seen.add(x); fr.append(x)
    return seen
def orbits(gens,n):
    par=list(range(n))
    def f(a):
        while par[a]!=a: par[a]=par[par[a]]; a=par[a]
        return a
    for g in gens:
        for i in range(n): par[f(i)]=f(g[i])
    O={}
    for i in range(n): O.setdefault(f(i),[]).append(i)
    return list(O.values())
for (N,t,k,v) in json.load(open(sys.argv[1])):
    q={7:2,13:3,31:5}[N]; f=v-N
    n,_,Hg,Ks,_=groups(q,f)
    from psl3 import pg2,norm
    pts=pg2(q); ix={p:i for i,p in enumerate(pts)}
    w=next(g for g in range(1,q) if len({pow(g,e,q) for e in range(q-1)})==q-1)
    def mat(M): return tuple([ix[norm(tuple(sum(M[r][c]*x[c] for c in range(3))%q for r in range(3)),q)] for x in pts]+list(range(n,v)))
    P3=[[0,1,0],[0,0,1],[1,0,0]]; P2=[[0,1,0],[1,0,0],[0,0,1]]; D=[[w,0,0],[0,1,0],[0,0,1]]
    U=[mat(P3),mat(P2),mat(D)]
    O=orbits(U,v); sizes=[len(o) for o in O]
    # subset sum
    reach={0:[]}
    for i,s in enumerate(sizes):
        for a,L in list(reach.items()):
            if a+s<=k and a+s not in reach: reach[a+s]=L+[i]
    ok=k in reach
    print(N,t,k,v,'certified' if ok else 'no', sorted(sizes))
