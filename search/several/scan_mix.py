"""Theorem 4.1 on X=F_p (+inf), H=<sigma,mu_d>, base designs with stabilisers A_e=<sigma,mu^(d/e)> (index d/e)
or B_e=<mu^(d/e)> (index pd/e); indices sum to N, at most 3 base designs, at least one of B type (else covered by Thm 3)."""
import sys, csv, time, itertools
from math import comb
from genls import solve, verify
def primroot(p):
    for g in range(2,p):
        if len({pow(g,i,p) for i in range(p-1)})==p-1: return g
def isprime(q): return q>1 and all(q%d for d in range(2,int(q**.5)+1))
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a'); MAXM=int(sys.argv[4]) if len(sys.argv)>4 else 3
jobs=[]
for row in csv.DictReader(open('../../table/existence_table.csv')):
    if row['status']!='open': continue
    N,t,k,v=map(int,(row['N'],row['t'],row['k'],row['v']))
    for inf in (0,1):
        p=v-inf
        if not isprime(p) or p<5: continue
        for d in range(1,p):
            if (p-1)%d: continue
            if comb(v,k)//(p*d)>CAP: continue
            subs=[]
            for e in range(1,d+1):
                if d%e==0: subs.append(('A',e,d//e)); subs.append(('B',e,p*d//e))
            for m in range(1,MAXM+1):
                for combo in itertools.combinations_with_replacement(subs,m):
                    if sum(c[2] for c in combo)==N and (m>1 or combo[0][0]=='B'):
                        jobs.append((comb(v,k)//(p*d),N,t,k,v,p,inf,d,combo))
jobs.sort(); print(len(jobs),"jobs",flush=True)
done=set()  # PREFILTER + RESUME
import os,re as _re
from filters import n1_n4
from necessary import count_sub
def _primes(n): return [q for q in range(2,n+1) if n%q==0 and all(q%d for d in range(2,q))]
def n2_closed(N,t,k,v,p,inf,d,combo,lam):
    for q in _primes(d):
        phi=sum(c[2] for c in combo if c[0]=='A' and c[1]%q==0)
        if any(c[0]=='B' for c in combo): return True
        ones=1+inf; qs=(p-1)//q
        for j in range(t//q+1):
            i=t-j*q
            if i>ones or j>qs: continue
            cnt=count_sub([1]*(ones-i)+[q]*(qs-j),k-t)
            if (cnt-phi*lam)%q or cnt<phi*(lam%q): return False
    return True
prev=set()
if os.path.exists(sys.argv[3]):
    for l in open(sys.argv[3]):
        prev.add(l.split(' SAT ok')[0].split(' UNSAT')[0].split(' UNKNOWN')[0].split(' INFEASIBLE')[0].split(' FILTER')[0].strip())
        if 'SAT ok= True' in l: done.add(tuple(map(int,_re.match(r'LS\[(\d+)\]\((\d+),(\d+),(\d+)\)',l).groups())))
print(len(prev),'already tried',flush=True)
for _,N,t,k,v,p,inf,d,combo in jobs:
    key=(N,t,k,v)
    if key in done: continue
    g=primroot(p); a=pow(g,(p-1)//d,p)
    def aff(m,b): return tuple([(m*x+b)%p for x in range(p)]+([p] if inf else []))
    sig=aff(1,1); H=[sig]+([aff(a,0)] if d>1 else [])
    Ks=[]
    for typ,e,_ in combo:
        mu=[aff(pow(a,d//e,p),0)] if e>1 else []
        Ks.append(([sig] if typ=='A' else [])+mu)
    lam=comb(v-t,k-t)//N; t0=time.time()
    tag=f"LS[{N}]({t},{k},{v}) X=F{p}{'+inf' if inf else ''} d={d} combo={combo}"
    if tag in prev: continue
    if not n2_closed(N,t,k,v,p,inf,d,combo,lam): print(tag,"FILTER N2",file=out,flush=True); continue
    f=n1_n4(v,t,k,lam,H,Ks,[c[2] for c in combo],tl=20)
    if f in ("N1","N4"): print(tag,"FILTER",f,file=out,flush=True); continue
    try: st,sol,info=solve(v,t,k,lam,H,Ks,tl=TL)
    except Exception as ex: print("ERR",key,ex,file=out,flush=True); continue
    tag=f"LS[{N}]({t},{k},{v}) X=F{p}{'+inf' if inf else ''} d={d} combo={combo}"
    if st=="SAT":
        ok=verify(v,t,k,lam,H,Ks,sol); print(tag,"SAT ok=",ok,"H=",H,"Ks=",Ks,"sol=",sol,file=out,flush=True)
        if ok: done.add(key)
    else: print(tag,st,f"[{time.time()-t0:.0f}s]",file=out,flush=True)
