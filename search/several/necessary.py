"""Necessary conditions N1 (stabiliser) and N2 (mod-q fixed-structure) for Theorem 4.1 data (H; K_1..K_m)."""
from itertools import combinations
from math import comb
from fix1 import gen_group
def cycles(g):
    seen=set(); cs=[]
    for x in range(len(g)):
        if x in seen: continue
        c=[x]; seen.add(x); y=g[x]
        while y!=x: c.append(y); seen.add(y); y=g[y]
        cs.append(len(c))
    return cs
def inv(g):
    r=[0]*len(g)
    for i,x in enumerate(g): r[x]=i
    return tuple(r)
def comp(a,b): return tuple(b[a[i]] for i in range(len(a)))  # a then b
def isprime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
def order(g):
    from math import lcm
    o=1
    for c in cycles(g): o=lcm(o,c)
    return o
def count_sub(lens, target):
    # number of sub-multisets (choice of distinct cycles) with total length target
    dp=[0]*(target+1); dp[0]=1
    for L in lens:
        for s in range(target,L-1,-1): dp[s]+=dp[s-L]
    return dp[target]
def N2(v,t,k,lam,Hel,Ksets):
    """returns (True,None) or (False, witness)"""
    seen=set()
    for g in Hel:
        q=order(g)
        if not isprime(q): continue
        ct=tuple(sorted(cycles(g)))
        # phi(g): number of designs fixed = sum_i #{h in H: h g h^-1 in K_i}/|K_i|
        phi=0
        for K in Ksets:
            c=sum(1 for h in Hel if comp(comp(inv(h),g),h) in K)
            phi+=c//len(K)
        key=(ct,phi)
        if key in seen: continue
        seen.add(key)
        ones=ct.count(1); qs=ct.count(q)
        # g-invariant t-sets: i fixed points + j q-cycles, i+jq=t
        for j in range(t//q+1):
            i=t-j*q
            if i>ones or j>qs: continue
            lens=[1]*(ones-i)+[q]*(qs-j)
            cnt=count_sub(lens,k-t)
            if (cnt-phi*lam)%q or cnt<phi*(lam%q): return False,dict(cycletype=ct,phi=phi,T=(i,j),count=cnt)
    return True,None
def check(v,t,k,lam,Hgens,Ks):
    Hel=gen_group(Hgens,v) if Hgens else {tuple(range(v))}
    Ksets=[frozenset(gen_group(K,v)) if K else frozenset({tuple(range(v))}) for K in Ks]
    return N2(v,t,k,lam,Hel,Ksets)
