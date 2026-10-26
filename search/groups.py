from core import *
from gf import GF
def ASL(q, m=1):
    """x->a x^(p^(j n/m)) + b, a nonzero square; returns generators as perms of range(q)."""
    F=GF(q); gens=[]
    g=next(a for a in range(2,q) if len({1}|set(_pw(F,a)))==q-1) if q>3 else 2
    sq=F.mulT[g][g]
    gens.append(tuple(F.add[x][1] for x in range(q)))
    if F.n>1:
        for i in range(1,F.n): gens.append(tuple(F.add[x][F.p**i] for x in range(q)))
    gens.append(tuple(F.mulT[sq][x] for x in range(q)))
    if m>1:
        e=F.n//m; 
        def fr(x):
            for _ in range(e): x=F.frob(x)
            return x
        gens.append(tuple(fr(x) for x in range(q)))
    return gens
def _pw(F,a):
    y=a; out=[a]
    while y!=1: y=F.mulT[y][a]; out.append(y)
    return out
def PSL2(q, m=1, pgl=False):
    """acting on range(q+1), point q = infinity."""
    F=GF(q); INF=q
    def mob(a,b,c,d):
        out=[]
        for x in range(q):
            num=F.add[F.mulT[a][x]][b]; den=F.add[F.mulT[c][x]][d]
            out.append(INF if den==0 else F.mulT[num][F.inv[den]])
        out.append(INF if c==0 else F.mulT[a][F.inv[c]])
        return tuple(out)
    g=next(a for a in range(2,q) if len(set(_pw(F,a)))==q-1) if q>3 else 2
    gens=[mob(1,1,0,1), mob(0,F.neg[1],1,0), mob(g if pgl else F.mulT[g][g],0,0,1)]
    if F.n>1:
        for i in range(1,F.n): gens.append(mob(1,F.p**i,0,1))
    if m>1:
        e=F.n//m
        def fr(x):
            if x==INF: return INF
            for _ in range(e): x=F.frob(x)
            return x
        gens.append(tuple(fr(x) for x in range(q+1)))
    return gens
def AGammaL1(q, full=True):
    F=GF(q)
    g=next(a for a in range(2,q) if len(set(_pw(F,a)))==q-1)
    gens=[tuple(F.add[x][F.p**i] for x in range(q)) for i in range(F.n)]
    gens.append(tuple(F.mulT[g][x] for x in range(q)))
    if full: gens.append(tuple(F.frob(x) for x in range(q)))
    return gens
