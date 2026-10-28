"""Theorem 4.1 on X = PG(1,q), v=q+1, H in {PGL(2,q), PSL(2,q) [q=3 mod 4], PGammaL(2,q) [q=p^m,m>1]}.
Stabilisers of small index: H, intermediate groups PSL<=L<=PGammaL, Borel-type F_q x| C_e (x| Frobenius)."""
import sys, csv, time, itertools
from math import comb, gcd
from genls2 import solve, verify
def ppow(q):
    for p in range(2,q+1):
        if q%p==0:
            m=0;r=q
            while r%p==0: r//=p; m+=1
            return (p,m) if r==1 else None
from scan_frob_gf import gf
def build(q):
    p,m=ppow(q); mul,add,w=gf(p,m) if m>1 else (lambda a,b:a*b%q, lambda a,b:(a+b)%q, next(g for g in range(2,q) if len({pow(g,i,q) for i in range(q-1)})==q-1))
    INF=q
    def inv(a):
        for b in range(1,q):
            if mul(a,b)==1: return b
    def pw(a,e):
        r=1
        for _ in range(e): r=mul(r,a)
        return r
    neg=lambda a: next(b for b in range(q) if add(a,b)==0)
    def perm(f): return tuple(f(x) for x in range(q+1))
    tr=lambda b: perm(lambda x: INF if x==INF else add(x,b))
    mu=lambda c: perm(lambda x: INF if x==INF else mul(c,x))
    m1=neg(1)
    iota=perm(lambda x: 0 if x==INF else (INF if x==0 else mul(m1,inv(x))))   # x -> -1/x (in PSL)
    iota2=perm(lambda x: 0 if x==INF else (INF if x==0 else inv(x)))          # x -> 1/x
    frob=perm(lambda x: INF if x==INF else pw(x,p))
    T=[tr(p**i) for i in range(m)]
    G={}
    n=(q-1)
    G['PGL']=(T+[mu(w),iota2], q*(q*q-1))
    G['PSL']=(T+[mu(pw(w,2)),iota], q*(q*q-1)//gcd(2,q-1))
    if m>1:
        G['PGamL']=(T+[mu(w),iota2,frob], m*q*(q*q-1))
        G['PSigL']=(T+[mu(pw(w,2)),iota,frob], m*q*(q*q-1)//gcd(2,q-1))
    for e in range(1,n+1):
        if n%e==0: G[f'B{e}']=(T+([mu(pw(w,n//e))] if e>1 else []), q*e)
        if n%e==0 and m>1: G[f'BF{e}']=(T+([mu(pw(w,n//e))] if e>1 else [])+[frob], q*e*m)
    return G
