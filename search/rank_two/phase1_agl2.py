"""Phase 1: for every target (N,t,k,v) and every group H of families.py on v points, list the pairs of
subgroup classes (K1,K2) with [H:K1]+[H:K2]=N that pass the derangement and congruence conditions."""
import sys, json, time
from math import comb
from pgrp import Group
from agl2 import agl2 as all_groups
from rank2 import cycle_type, coset_fix, cheap
targets=json.load(open(sys.argv[1])); out=open(sys.argv[2],'w'); MAXO=int(sys.argv[3]) if len(sys.argv)>3 else 6000
byv={}
for x in targets: byv.setdefault(x[-1],[]).append(tuple(x[-4:]))
for v in sorted(byv):
    for name,gens in all_groups(v,MAXO):
        t0=time.time()
        G=Group(gens,v)
        Ns={T[0] for T in byv[v]}
        Ns={N for N in Ns if G.n%N==0 or any(G.n%a==0 and G.n%(N-a)==0 for a in range(1,N))}
        if not Ns: continue
        try: C=G.subgroup_classes(max_classes=4000)
        except RuntimeError: print("SKIP",v,name,file=sys.stderr,flush=True); continue
        byidx={}
        for K in C:
            byidx.setdefault(G.n//len(K),[]).append(K)
        pairs={}
        for N in Ns:
            L=[(K,) for K in byidx.get(N,[])]
            for a in sorted(byidx):
                if a<N-a and (N-a) in byidx:
                    for K1 in byidx[a]:
                        for K2 in byidx[N-a]: L.append((K1,K2))
            if L: pairs[N]=L
        if not pairs: continue
        cts=[cycle_type(G.perm(i)) for i in range(G.n)]
        phic={}
        for N,L in pairs.items():
            for Ks in L:
                for K in Ks:
                    if K not in phic: phic[K]=coset_fix(G,K)
                phi=sum(phic[K] for K in Ks)
                for T in byv[v]:
                    if T[0]!=N: continue
                    _,t,k,_=T; lam=comb(v-t,k-t)//N
                    if cheap(G,phi,t,k,v,lam,cts)=="ok":
                        rec=dict(target=T,H=name,Hgens=[G.perm(i) for i in G.gens_idx],
                                 Kgens=[[G.perm(i) for i in G.gens_of(K)] for K in Ks],
                                 idx=[G.n//len(K) for K in Ks],Kord=[len(K) for K in Ks],Horder=G.n)
                        print(json.dumps(rec),file=out,flush=True)
        print(v,name,G.n,len(C),f"{time.time()-t0:.1f}s",file=sys.stderr,flush=True)
