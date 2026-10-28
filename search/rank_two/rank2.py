"""Rank-two (and rank-three) frame search: H from families.py, K_i over all subgroup classes of H."""
import sys, json, time, os
from math import comb
from itertools import combinations_with_replacement
from pgrp import Group
from families import all_groups
from filters import n1_n4
from genls2 import solve, verify
GCACHE={}
def group(gens,v):
    key=(v,tuple(gens))
    if key not in GCACHE:
        G=Group(gens,v); G.classes=None; GCACHE[key]=G
    return GCACHE[key]
def cycle_type(perm):
    seen=[False]*len(perm); L=[]
    for x in range(len(perm)):
        if not seen[x]:
            l=0;y=x
            while not seen[y]: seen[y]=True; y=perm[y]; l+=1
            L.append(l)
    return L
def fixes_kset(L,k):
    s={0}
    for l in L: s|={a+l for a in s if a+l<=k}
    return k in s
def coset_fix(G,K):
    """phi_K(g) for all g: number of right cosets Kh fixed by g."""
    Ks=K; reps=[]; cov=set()
    for h in range(G.n):
        if h in cov: continue
        reps.append(h); cov|={int(G.M[k,h]) for k in Ks}
    import numpy as np
    phi=np.zeros(G.n,dtype=np.int32)
    Kset=set(Ks)
    for h in reps:
        hi=G.inv[h]
        # g fixes Kh iff h g h^-1 in K iff g in h^-1 K h
        for k in Ks: phi[int(G.M[G.M[hi,k],h])]+=1
    return phi
def cheap(G,phi,t,k,v,lam,cts):
    for g in range(1,G.n):
        if phi[g]==0 and fixes_kset(cts[g],k): return "DER"
    for g in range(1,G.n):
        L=cts[g]; q=max(L)
        if q==1 or any(l not in (1,q) for l in L): continue
        if any(q%d==0 and d not in(1,q) for d in range(2,q)): continue  # prime order only
        f0=L.count(1); c=L.count(q); ph=int(phi[g]); r=lam%q
        for y in range(0,t//q+1):
            x=t-y*q
            if x>f0 or y>c: continue
            cnt=sum(comb(f0-x,x2)*comb(c-y,(k-t-x2)//q) for x2 in range(0,k-t+1) if (k-t-x2)%q==0 and (k-t-x2)//q>=0)
            if (cnt-ph*lam)%q or cnt<ph*r: return "CONG"
    return "ok"
def run(target,out,maxorder=3000,rank=2,tl=600,ftl=60,logf=None):
    N,t,k,v=target; lam=comb(v-t,k-t)//N
    for name,gens in all_groups(v,maxorder):
        tg=f"LS[{N}]({t},{k},{v}) H={name}"
        if tg in DONE: continue
        t0=time.time()
        G=group(gens,v)
        if G.classes is None:
            try: G.classes=G.subgroup_classes()
            except RuntimeError: G.classes=[]
        subs=[K for K in G.classes if G.n//len(K)<=N-1]
        byidx={}
        for K in subs: byidx.setdefault(G.n//len(K),[]).append(K)
        cts=[cycle_type(G.perm(i)) for i in range(G.n)]
        tried=0; passed=0; res=[]
        idxs=sorted(byidx)
        combos=[]
        if rank==2:
            for a in idxs:
                b=N-a
                if b in byidx and a<b:
                    for K1 in byidx[a]:
                        for K2 in byidx[b]: combos.append((K1,K2))
        for Ks in combos:
            tried+=1
            phi=sum(coset_fix(G,K) for K in Ks)
            c=cheap(G,phi,t,k,v,lam,cts)
            if c!="ok": continue
            passed+=1
            Kg=[[G.perm(i) for i in G.gens_of(K)] for K in Ks]
            Hg=[G.perm(i) for i in G.gens_idx]
            ids=[G.n//len(K) for K in Ks]
            fl=n1_n4(v,t,k,lam,Hg,Kg,ids,tl=ftl)
            if fl not in ("pass","N4?"): res.append(("FILTER",fl,ids,[len(K) for K in Ks])); continue
            st,sol,info=solve(v,t,k,lam,Hg,Kg,ids,tl=tl,workers=2)
            rec=dict(target=target,H=name,Hgens=Hg,Kgens=Kg,idx=ids,status=st,info=info)
            if st=="SAT":
                rec['ok']=verify(v,t,k,lam,Hg,Kg,sol,N); rec['sol']=sol
            print(json.dumps(rec),file=out,flush=True)
            res.append((st,ids,[len(K) for K in Ks]))
            if st=="SAT" and rec['ok']: break
        print(tg,f"|H|={G.n} classes={len(G.classes)} pairs={tried} cheap-pass={passed}",res,f"{time.time()-t0:.0f}s",file=logf,flush=True)
        if any(r[0]=="SAT" for r in res): return True
    return False
if __name__=="__main__":
    targets=json.load(open(sys.argv[1])); out=open(sys.argv[2],'a'); logf=open(sys.argv[3],'a')
    DONE=set(l.split(' |H|')[0] for l in open(sys.argv[3])) if os.path.exists(sys.argv[3]) else set()
    for T in targets:
        r=run(tuple(T[-4:]),out,logf=logf)
        print("TARGET",T,"FOUND" if r else "none",file=logf,flush=True)
