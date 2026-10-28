"""Dedicated m=1 frame: X = F_p u F (|F|=f=v-p < p), H = (T x| C_e) x Sym(F), K = C_e x Sym(F) (point stabiliser), [H:K]=p=N.
Rank-one, non-regular P = C_p x| C_e. Requires f < k mod p (shifter condition)."""
import sys, csv, time, os
from math import comb
from genls2 import solve, verify
from filters import n1_n4
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a')
def primroot(p): return next(g for g in range(2,p) if len({pow(g,i,p) for i in range(p-1)})==p-1)
def isprime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
rows=[r for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','table','existence_table.csv'))) if r['status']=='open']
T=sorted({(comb(int(r['v']),int(r['k'])),int(r['N']),int(r['t']),int(r['k']),int(r['v'])) for r in rows
   if isprime(int(r['N'])) and int(r['N'])>=5 and 0<=int(r['v'])-int(r['N'])<int(r['k'])%int(r['N']) and comb(int(r['v']),int(r['k']))<=CAP})
print(len(T),"targets",flush=True)
prev=set(l.split(' ::')[0] for l in open(sys.argv[3])) if os.path.exists(sys.argv[3]) else set()
for _,p,t,k,v in T:
    f=v-p; lam=comb(v-t,k-t)//p; g=primroot(p)
    P=lambda fn: tuple(fn(x) for x in range(v))
    sig=P(lambda x:(x+1)%p if x<p else x)
    Fg=[]
    if f>=2: Fg.append(P(lambda x: x if x<p else p+((x-p+1)%f)))
    if f>=3: Fg.append(P(lambda x: {p:p+1,p+1:p}.get(x,x)))
    for e in sorted([e for e in range(2,p) if (p-1)%e==0],reverse=True)[:3]:
        a=pow(g,(p-1)//e,p); mu=P(lambda x: a*x%p if x<p else x)
        H=[sig,mu]+Fg; K=[mu]+Fg
        tag=f"LS[{p}]({t},{k},{v}) F{p}+{f} e={e}"
        if tag in prev: continue
        try:
            fl=n1_n4(v,t,k,lam,H,[K],[p],tl=20)
            if fl not in ("pass","N4?"): print(tag,":: FILTER",fl,file=out,flush=True); continue
            t0=time.time(); st,sol,info=solve(v,t,k,lam,H,[K],[p],tl=TL,workers=1)
        except ValueError as ex: print(tag,":: ERR",ex,file=out,flush=True); continue
        if st=="SAT":
            ok=verify(v,t,k,lam,H,[K],sol,p); print(tag,":: SAT ok=",ok,"H=",H,"Ks=",[K],"sol=",sol,file=out,flush=True)
            if ok: break
        else: print(tag,":: ",st,info,f"[{time.time()-t0:.0f}s]",file=out,flush=True)
