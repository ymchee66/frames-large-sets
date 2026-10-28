"""m=1 search (one base design, P transitive of prime degree N) over three families of frames:
 A/B: X=F_q (+inf), H=<T',C_{eN},psi>, K=<T',C_e,psi>, psi a Frobenius power or none, T' in {T, translations of a subfield F_{p^j} containing C_{eN}};
 C:   X=PG(1,q), H=<L,phi>, K=<L,phi^N>, L in {PSL(2,q),PGL(2,q)}, N | a.
Necessary condition (Theorem 3.1 of the note): an element of order N mapping to an N-cycle has v mod N < k mod N fixed points."""
import sys, csv, time, os
from math import comb
from genls2 import solve, verify
from filters import n1_n4
from scan_frob_gf import gf
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a')
def ppow(q):
    for p in range(2,q+1):
        if q%p==0:
            m=0;r=q
            while r%p==0: r//=p; m+=1
            return (p,m) if r==1 else None
def isprime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
FIELD={}
def field(q):
    if q not in FIELD:
        p,a=ppow(q)
        if a==1:
            g=next(g for g in range(2,q) if len({pow(g,i,q) for i in range(q-1)})==q-1)
            FIELD[q]=(p,a,lambda x,y:x*y%q,lambda x,y:(x+y)%q,g)
        else:
            mul,add,w=gf(p,a); FIELD[q]=(p,a,mul,add,w)
    return FIELD[q]
def affine_frames(q,inf,N):
    p,a,mul,add,w=field(q)
    def pw(x,e):
        r=1
        for _ in range(e): r=mul(r,x)
        return r
    ext=[q] if inf else []
    perm=lambda f: tuple([f(x) for x in range(q)]+ext)
    tr=lambda b: perm(lambda x: add(x,b)); mu=lambda c: perm(lambda x: mul(c,x)); frob=lambda c: perm(lambda x: pw(x,p**c))
    res=[]
    for eN in range(N,q,N):
        if (q-1)%eN: continue
        e=eN//N
        cN=pw(w,(q-1)//eN); ce=pw(cN,N)
        for c in [None]+[c for c in range(1,a) if a%c==0]:     # psi = phi^c of order a/c
            psi=[frob(c)] if c else []
            # translation part: full T, or subfield F_{p^j} with C_eN inside F_{p^j}^* and psi preserving it
            for j in [a]+[j for j in range(1,a) if a%j==0 and (p**j-1)%eN==0]:
                if j==a: Tp=[tr(p**i) for i in range(a)]
                else:
                    z=pw(w,(q-1)//(p**j-1)); Tp=[tr(1)]+[tr(pw(z,i)) for i in range(1,j)]
                H=Tp+[mu(cN)]+psi; K=Tp+([mu(ce)] if e>1 else [])+psi
                res.append((f"F{q}{'+inf' if inf else ''} T'=F_{p}^{j} e={e} psi={'phi^%d'%c if c else '-'}",H,K))
    return res
def proj_frames(q,N):
    from psl_build import build
    p,a=ppow(q)
    if a%N: return []
    G=build(q); res=[]
    phi=G['PGamL'][0][-1]
    def power(g,n):
        r=tuple(range(len(g)))
        for _ in range(n): r=tuple(g[r[i]] for i in range(len(g)))
        return r
    for L in ('PSL','PGL'):
        Lg=G[L][0][:-1] if False else [g for g in G[L][0]]
        res.append((f"PG(1,{q}) H={L}.<phi> K={L}.<phi^{N}>", Lg+[phi], Lg+[power(phi,N)]))
    return res

rows=[r for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','table','existence_table.csv'))) if r['status'] in ('open','exists')]
jobs=[]
for r in rows:
    N,t,k,v=map(int,(r['N'],r['t'],r['k'],r['v']))
    if N<3 or comb(v,k)>CAP: continue
    for inf in (0,1):
        q=v-inf; pm=ppow(q)
        if pm and pm[1]>1 and (q-1)%N==0 and 1+inf<k%N: jobs.append((comb(v,k),N,t,k,v,'B',q,inf,r['status']))
    q=v-1; pm=ppow(q)
    if pm and pm[1]%N==0 and pm[0]**(pm[1]//N)+1<k%N: jobs.append((comb(v,k),N,t,k,v,'C',q,0,r['status']))
jobs.sort(); print(len(jobs),"jobs",sorted({j[1:5]+(j[5],j[8]) for j in jobs}),flush=True)
prev=set(l.split(' ::')[0] for l in open(sys.argv[3])) if os.path.exists(sys.argv[3]) else set()
for _,N,t,k,v,fam,q,inf,stt in jobs:
    lam=comb(v-t,k-t)//N
    frames=[f for f in affine_frames(q,inf,N) if 'psi=phi' in f[0]] if fam=='B' else proj_frames(q,N)
    for desc,H,K in frames:
        tag=f"LS[{N}]({t},{k},{v}) [{stt}] {desc}"
        if tag in prev: continue
        try:
            f=n1_n4(v,t,k,lam,H,[K],[N],tl=30)
        except Exception as ex:
            print(tag,":: ERR",ex,file=out,flush=True); continue
        if f not in ("pass","N4?"):
            print(tag,":: FILTER",f,file=out,flush=True); continue
        t0=time.time()
        try: st,sol,info=solve(v,t,k,lam,H,[K],[N],tl=TL)
        except ValueError as ex:
            print(tag,":: ERR",ex,file=out,flush=True); continue
        if st=="SAT":
            ok=verify(v,t,k,lam,H,[K],sol,N)
            print(tag,":: SAT ok=",ok,"H=",H,"Ks=",[K],"sol=",sol,file=out,flush=True)
        else: print(tag,":: ",st,info,f"[{time.time()-t0:.0f}s]",file=out,flush=True)
