"""Candidate groups H on v points: affine (F_q + F), projective (PG(1,q) + F), cyclic-multiplier (Z_n + F).
F-part: fixed pointwise, full Sym(F) as direct factor, or a permutation of F coupled to the multiplier."""
from math import gcd
from scan_frob_gf import gf
def ppow(q):
    for p in range(2,q+1):
        if q%p==0:
            m=0;r=q
            while r%p==0: r//=p; m+=1
            return (p,m) if r==1 else None
def field(q):
    p,a=ppow(q)
    if a==1:
        g=next(g for g in range(2,q) if len({pow(g,i,q) for i in range(q-1)})==q-1) if q>2 else 1
        return p,a,(lambda x,y:x*y%q),(lambda x,y:(x+y)%q),g
    mul,add,w=gf(p,a); return p,a,mul,add,w
def fparts(base,v):
    """F = points base..v-1"""
    f=v-base
    out=[("F1",[],None)]
    if f>=2:
        sw=tuple(x if x<base else ({base:base+1,base+1:base}.get(x,x)) for x in range(v))
        cy=tuple(x if x<base else base+((x-base+1)%f) for x in range(v))
        out.append(("FSym",[cy]+([sw] if f>=3 else []),None))
        out.append(("Fsw",[],sw))
        if f>=3: out.append(("Fcy",[],cy))
    return out
def comp(a,b): return tuple(b[a[i]] for i in range(len(a)))
def affine(v,maxorder):
    res=[]
    for f in range(0,4):
        q=v-f
        if q<3 or not ppow(q): continue
        p,a,mul,add,w=field(q)
        def pw(x,e):
            r=1
            for _ in range(e): r=mul(r,x)
            return r
        perm=lambda fn: tuple([fn(x) for x in range(q)]+list(range(q,v)))
        T=[perm(lambda x,b=p**i: add(x,b)) for i in range(a)]
        for e in [d for d in range(1,q) if (q-1)%d==0]:
            mu=perm(lambda x,c=pw(w,(q-1)//e): mul(c,x))
            for c in [None]+[c for c in range(1,a) if a%c==0]:
                psi=[perm(lambda x,c=c: pw(x,p**c))] if c else []
                for fn,Sg,pi in fparts(q,v):
                    if pi is not None and e==1: continue
                    m=comp(mu,pi) if pi is not None else mu
                    gens=T+([m] if e>1 else [])+psi+Sg
                    order=q*e*((a//c) if c else 1)*({"F1":1,"FSym":1,"Fsw":1,"Fcy":1}[fn])
                    if fn=="FSym": order*= {2:2,3:6}[v-q]
                    if order>maxorder: continue
                    res.append((f"AF{q}+{f} e={e}{' phi^%d'%c if c else ''} {fn}",gens))
    return res
def projective(v,maxorder):
    from psl_build import build
    res=[]
    for f in range(0,4):
        q=v-1-f
        if q<4 or not ppow(q): continue
        G=build(q)
        for name in ['PSL','PGL','PSigL','PGamL']:
            if name not in G: continue
            gens,order=G[name]
            if order>maxorder: continue
            gens=[tuple(list(g)+list(range(q+1,v))) for g in gens]
            for fn,Sg,pi in fparts(q+1,v):
                if pi is not None: continue      # coupling needs a quotient map; skip
                o=order*({2:2,3:6}.get(f,1) if fn=="FSym" else 1)
                if o>maxorder: continue
                res.append((f"P{name}({q})+{f} {fn}",gens+Sg))
    return res
def cyclic(v,maxorder):
    res=[]
    for f in range(0,4):
        n=v-f
        if n<4 or ppow(n): continue   # composite non prime power
        units=[u for u in range(1,n) if gcd(u,n)==1]
        def ordr(u):
            o=1;x=u
            while x!=1: x=x*u%n; o+=1
            return o
        seen=set()
        for u in units:
            sub=frozenset(pow(u,i,n) for i in range(ordr(u)))
            if sub in seen: continue
            seen.add(sub)
            perm=lambda fn: tuple([fn(x) for x in range(n)]+list(range(n,v)))
            sig=perm(lambda x:(x+1)%n); mu=perm(lambda x,u=u:u*x%n)
            for fn,Sg,pi in fparts(n,v):
                if pi is not None and u==1: continue
                m=comp(mu,pi) if pi is not None else mu
                gens=[sig]+([m] if u!=1 else [])+Sg
                o=n*len(sub)*({2:2,3:6}.get(f,1) if fn=="FSym" else 1)
                if o>maxorder: continue
                res.append((f"Z{n}+{f} u={u}(ord {len(sub)}) {fn}",gens))
    return res
def all_groups(v,maxorder=3000):
    return affine(v,maxorder)+projective(v,maxorder)+cyclic(v,maxorder)
if __name__=="__main__":
    import sys
    for nm,g in all_groups(int(sys.argv[1])): print(nm,len(g))

# ---- higher-dimensional affine groups T x| L, L running over subgroup classes of GL(n,p) or GL(2,4)
from itertools import product as _prod
def _affine_n(n,p,v,maxorder):
    from pgrp import Group
    q=p**n
    vecs=list(_prod(range(p),repeat=n)); ix={x:i for i,x in enumerate(vecs)}
    def lin(M): return tuple(ix[tuple(sum(M[r][c]*x[c] for c in range(n))%p for r in range(n))] for x in vecs)
    def trans(b): return tuple(ix[tuple((x[i]+b[i])%p for i in range(n))] for x in vecs)
    # generators of GL(n,p): elementary matrices and a diagonal generator
    gens=[]
    for i in range(n):
        for j in range(n):
            if i!=j:
                M=[[int(r==c) for c in range(n)] for r in range(n)]; M[i][j]=1; gens.append(lin(M))
    g=next(g for g in range(1,p) if p==2 or len({pow(g,e,p) for e in range(p-1)})==p-1)
    M=[[int(r==c) for c in range(n)] for r in range(n)]; M[0][0]=g; gens.append(lin(M))
    GL=Group(gens,q)
    Ls=GL.subgroup_classes(max_classes=3000)
    T=[trans(tuple(int(i==j) for j in range(n))) for i in range(n)]
    ext=lambda pm: tuple(list(pm)+list(range(q,v)))
    res=[]
    for L in Ls:
        Lg=[GL.perm(i) for i in GL.gens_of(L)]
        for fn,Sg,pi in fparts(q,v):
            if pi is not None: continue
            o=q*len(L)*({2:2,3:6}.get(v-q,1) if fn=="FSym" else 1)
            if o>maxorder: continue
            res.append((f"AGL({n},{p})+{v-q} |L|={len(L)} {fn}",[ext(x) for x in T]+[ext(x) for x in Lg]+Sg))
    return res
def affine_n(v,maxorder):
    res=[]
    for (n,p) in [(3,2),(2,3)]:
        q=p**n
        if 0<=v-q<=3: res+=_affine_n(n,p,v,maxorder)
    return res
_old_all=all_groups
def all_groups(v,maxorder=3000):
    return _old_all(v,maxorder)+affine_n(v,maxorder)
