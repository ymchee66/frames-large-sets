"""Groups T x| L on F_q^2 (+ f fixed points), L over subgroup classes of GammaL(2,4) or GL(2,5)."""
from itertools import product
from pgrp import Group
from scan_frob_gf import gf
def field(q):
    if q==5: return (lambda a,b:a*b%5),(lambda a,b:(a+b)%5),2,5,1
    mul,add,w=gf(2,2); return mul,add,w,2,2
def agl2(v,maxorder=6000,fparts=None):
    res=[]
    for q in (4,5):
        n=q*q; f=v-n
        if f<0 or f>4: continue
        mul,add,w,p,e=field(q)
        vecs=list(product(range(q),repeat=2)); ix={x:i for i,x in enumerate(vecs)}
        def lin(M): return tuple(ix[(add(mul(M[0][0],x[0]),mul(M[0][1],x[1])),add(mul(M[1][0],x[0]),mul(M[1][1],x[1])))] for x in vecs)
        gens=[lin([[1,1],[0,1]]),lin([[1,0],[1,1]]),lin([[w,0],[0,1]])]
        if q==4: gens.append(tuple(ix[(mul(x[0],x[0]),mul(x[1],x[1]))] for x in vecs))   # Frobenius
        GL=Group(gens,n)
        Ls=GL.subgroup_classes(max_classes=4000)
        T=[tuple(ix[(add(x[0],b[0]),add(x[1],b[1]))] for x in vecs) for b in [(1,0),(0,1)]+([(w,0),(0,w)] if q==4 else [])]
        ext=lambda pm: tuple(list(pm)+list(range(n,v)))
        sw=tuple(list(range(n))+[n+1,n]+list(range(n+2,v))) if f>=2 else None
        cy=tuple(list(range(n))+[n+((i+1)%f) for i in range(f)]) if f>=2 else None
        for L in Ls:
            Lg=[GL.perm(i) for i in GL.gens_of(L)]
            base=[ext(x) for x in T]+[ext(x) for x in Lg]
            o=n*len(L)
            if o<=maxorder: res.append((f"AGL2({q})+{f} |L|={len(L)} F1",base))
            if f>=2:
                Sg=[cy]+([sw] if f>=3 else [])
                fs={2:2,3:6,4:24}[f]
                if o*fs<=maxorder: res.append((f"AGL2({q})+{f} |L|={len(L)} FSym",base+Sg))
    return res
if __name__=="__main__":
    import sys
    L=agl2(int(sys.argv[1])); print(len(L)); print([x[0] for x in L][:10])
