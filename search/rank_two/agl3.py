"""Groups T x| L on F_3^3 (+ f extra points), L over subgroup classes of GL(3,3) with 27|L| <= maxorder."""
from itertools import product
from pgrp import Group
import pickle, os
CACHE='/home/claude/ls/rank2/gl33_classes.pkl'
def gl33_classes(maxL):
    if os.path.exists(CACHE):
        d=pickle.load(open(CACHE,'rb'))
        if d['maxL']>=maxL: return d['vecs'],[L for L in d['Ls'] if len(L[0]) and len(L[1])<=maxL] if False else [x for x in d['Ls'] if x[1]<=maxL]
    p=3; n=3
    vecs=list(product(range(p),repeat=n)); ix={x:i for i,x in enumerate(vecs)}
    def lin(M): return tuple(ix[tuple(sum(M[r][c]*x[c] for c in range(n))%p for r in range(n))] for x in vecs)
    gens=[]
    for i in range(n):
        for j in range(n):
            if i!=j:
                M=[[int(r==c) for c in range(n)] for r in range(n)]; M[i][j]=1; gens.append(lin(M))
    M=[[int(r==c) for c in range(n)] for r in range(n)]; M[0][0]=2; gens.append(lin(M))
    G=Group(gens,27)
    print("GL(3,3) order",G.n,flush=True)
    Ls=G.subgroup_classes(max_classes=20000,max_order=maxL)
    out=[([G.perm(i) for i in G.gens_of(L)],len(L)) for L in Ls]
    pickle.dump(dict(maxL=maxL,vecs=vecs,Ls=out),open(CACHE,'wb'))
    return vecs,out
def agl3(v,maxorder=6000):
    f=v-27
    if f<0 or f>4: return []
    vecs,Ls=gl33_classes(maxorder//27)
    ix={x:i for i,x in enumerate(vecs)}
    T=[tuple(ix[tuple((x[i]+b[i])%3 for i in range(3))] for x in vecs) for b in [(1,0,0),(0,1,0),(0,0,1)]]
    ext=lambda pm: tuple(list(pm)+list(range(27,v)))
    sw=tuple(list(range(27))+[28,27]+list(range(29,v))) if f>=2 else None
    cy=tuple(list(range(27))+[27+((i+1)%f) for i in range(f)]) if f>=2 else None
    res=[]
    for Lg,o in Ls:
        base=[ext(x) for x in T]+[ext(x) for x in Lg]
        if 27*o<=maxorder: res.append((f"AGL3(3)+{f} |L|={o} F1",base))
        if f>=2:
            fs={2:2,3:6,4:24}[f]
            if 27*o*fs<=maxorder: res.append((f"AGL3(3)+{f} |L|={o} FSym",base+[cy]+([sw] if f>=3 else [])))
    return res
if __name__=="__main__":
    import time; t=time.time()
    L=agl3(27); print(len(L), time.time()-t)
