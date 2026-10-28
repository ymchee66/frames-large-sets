"""Family E: N=p prime, X=F_p u F (|F|=f=v-p, f < k mod p), H <= AGL(1,p) x Sym(F) with image T x| C_r,
K = elements whose F_p-part fixes 0. Kernel part S in {1, Sym(F)}; multiplier mu_r combined with pi in {id, swap, f-cycle}."""
import sys, csv, time, os
from math import comb
from genls2 import solve, verify
from filters import n1_n4
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a')
def primroot(p): return next(g for g in range(2,p) if len({pow(g,i,p) for i in range(p-1)})==p-1) if p>2 else 1
def isprime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
rows=[r for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','table','existence_table.csv'))) if r['status'] in ('open','exists')]
ST={(int(r['N']),int(r['t']),int(r['k']),int(r['v'])):r['status'] for r in rows}
T=sorted({(comb(int(r['v']),int(r['k'])),int(r['N']),int(r['t']),int(r['k']),int(r['v'])) for r in rows
   if isprime(int(r['N'])) and int(r['N'])>=3 and 0<=int(r['v'])-int(r['N'])<int(r['k'])%int(r['N']) and comb(int(r['v']),int(r['k']))<=CAP})
print(len(T),"targets",flush=True)
prev=set(l.split(' ::')[0] for l in open(sys.argv[3])) if os.path.exists(sys.argv[3]) else set()
for _,p,t,k,v in T:
    f=v-p; lam=comb(v-t,k-t)//p; g=primroot(p)
    P=lambda fn: tuple(fn(x) for x in range(v))
    comp=lambda a,b: tuple(b[a[i]] for i in range(v))
    sig=P(lambda x:(x+1)%p if x<p else x)
    Fcyc=P(lambda x: x if x<p else p+((x-p+1)%f)) if f>=2 else None
    Fsw=P(lambda x: {p:p+1,p+1:p}.get(x,x)) if f>=2 else None
    Sym=[g_ for g_ in ([Fcyc] if f>=2 else [])+([Fsw] if f>=3 else [])]
    found=False
    for r in sorted([r for r in range(2,p) if (p-1)%r==0],reverse=True):
        a=pow(g,(p-1)//r,p); mu=P(lambda x: a*x%p if x<p else x)
        pis=[("id",None)]
        if f>=2 and r%2==0: pis.append(("swap",Fsw))
        if f>=3 and r%f==0: pis.append(("cyc",Fcyc))
        for sname,S in (("S=1",[]),("S=Sym",Sym)) if f>=2 else (("S=1",[]),):
            for pn,pi in pis:
                if S and pi is not None: continue      # pi absorbed in Sym(F)
                m=comp(mu,pi) if pi else mu
                H=[sig,m]+S; K=[m]+S
                tag=f"LS[{p}]({t},{k},{v}) [{ST[(p,t,k,v)]}] E F{p}+{f} r={r} {sname} pi={pn}"
                if tag in prev: continue
                try:
                    fl=n1_n4(v,t,k,lam,H,[K],[p],tl=20)
                    if fl not in ("pass","N4?"): print(tag,":: FILTER",fl,file=out,flush=True); continue
                    t0=time.time(); st,sol,info=solve(v,t,k,lam,H,[K],[p],tl=TL,workers=1)
                except ValueError as ex: print(tag,":: ERR",ex,file=out,flush=True); continue
                if st=="SAT":
                    ok=verify(v,t,k,lam,H,[K],sol,p); print(tag,":: SAT ok=",ok,"H=",H,"Ks=",[K],"sol=",sol,file=out,flush=True)
                else: print(tag,":: ",st,info,f"[{time.time()-t0:.0f}s]",file=out,flush=True)
