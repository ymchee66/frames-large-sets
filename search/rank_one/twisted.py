"""Twisted frame: X = F_p u {inf1, inf2}; H = <sigma, mu_g * swap> (mu_a swaps the infinities iff a is a non-square),
K = <mu_g * swap> (point stabiliser of 0), [H:K]=p. One base design (m=1), P = AGL(1,p) acting on designs as on F_p."""
import sys, time
from math import comb
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from genls2 import solve, verify
from filters import n1_n4
def primroot(p): return next(g for g in range(2,p) if len({pow(g,i,p) for i in range(p-1)})==p-1)
out=open('twisted.log','a')
for spec in sys.argv[2:]:
    p,k=map(int,spec.split(','))
    TL=float(sys.argv[1]); v=p+2; t=2; lam=comb(v-t,k-t)//p; g=primroot(p)
    P=lambda fn: tuple(fn(x) for x in range(v))
    sig=P(lambda x:(x+1)%p if x<p else x)
    tw=P(lambda x: g*x%p if x<p else (p+1 if x==p else p))     # mu_g (g non-square) swaps the infinities
    H=[sig,tw]; K=[tw]
    t0=time.time()
    fl=n1_n4(v,t,k,lam,H,[K],[p],tl=30)
    if fl not in ("pass","N4?"):
        print(f"LS[{p}](2,{k},{v}) twisted :: FILTER {fl}",file=out,flush=True); print(p,k,"filter",fl,flush=True); continue
    st,sol,info=solve(v,t,k,lam,H,[K],[p],tl=TL,workers=1)
    ok=verify(v,t,k,lam,H,[K],sol,p) if sol else None
    print(f"LS[{p}](2,{k},{v}) twisted :: {st} ok={ok} {info} [{time.time()-t0:.0f}s]"+(f" H= {H} Ks= {[K]} sol= {sol}" if sol else ""),file=out,flush=True)
    print(p,k,st,ok,f"[{time.time()-t0:.0f}s]",flush=True)
