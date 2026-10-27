"""Theorem 4.1 with Frobenius-type groups A x| M on X = A (+inf):
 (Z) A = Z_n (n composite), M = <a> <= Z_n^*;  stabilisers <sigma^s, mu^(d/e)>, index s*d/e (s|n, e|d)
 (F) A = (F_q,+), q = p^m, m>=2, M = <w^((q-1)/d)>; stabilisers B x| M_e with B in {0, subfields containing M_e, F_q}, index (q/|B|)*(d/e)
 at most 3 base designs, indices summing to N."""
import sys, csv, time, itertools
from math import comb, gcd
from genls import solve, verify
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a'); MAXM=3
def isprime(q): return q>1 and all(q%d for d in range(2,int(q**.5)+1))
def ppow(q):
    for p in range(2,q+1):
        if q%p==0:
            m=0; r=q
            while r%p==0: r//=p; m+=1
            return (p,m) if r==1 else None
# ---------- GF(p^m) ----------
def gf(p,m):
    import itertools as it
    for c in it.product(range(p),repeat=m):
        f=list(c)+[1]   # monic x^m + ...
        q=p**m
        def mul(a,b):
            A=[(a//p**i)%p for i in range(m)]; B=[(b//p**i)%p for i in range(m)]
            P=[0]*(2*m-1)
            for i in range(m):
                for j in range(m): P[i+j]=(P[i+j]+A[i]*B[j])%p
            for d in range(2*m-2,m-1,-1):
                co=P[d]
                if co:
                    for i in range(m+1): P[d-m+i]=(P[d-m+i]-co*f[i])%p
            return sum(P[i]*p**i for i in range(m))
        def add(a,b): return sum((((a//p**i)+(b//p**i))%p)*p**i for i in range(m))
        # primitive element test: x (=p) generates
        x=p if m>1 else 2; y=1; seen=set()
        for _ in range(q-1): y=mul(y,x); seen.add(y)
        if len(seen)==q-1 and 0 not in seen: return mul,add,x
    raise ValueError
jobs=[]
rows=[r for r in csv.DictReader(open('../../table/existence_table.csv')) if r['status']=='open']
fields={}
for row in rows:
    N,t,k,v=map(int,(row['N'],row['t'],row['k'],row['v']))
    for inf in (0,1):
        n=v-inf
        if n<6: continue
        pm=ppow(n)
        if pm and pm[1]>=2:   # field case
            p,m=pm; q=n
            subf=[p**r for r in range(1,m) if m%r==0]
            for d in range(2,q):
                if (q-1)%d or comb(v,k)//(q*d)>CAP: continue
                subs=[]
                for e in range(1,d+1):
                    if d%e: continue
                    subs.append(('B0',e,q*d//e)); subs.append(('BA',e,d//e))
                    for s in subf:
                        if ((q-1)//(d//e)*0==0) and (s-1)%e==0 and ( (q-1)//d * (d//e) ) % ((q-1)//(s-1))==0:
                            subs.append((f'BF{s}',e,(q//s)*(d//e)))
                for mm in range(1,MAXM+1):
                    for combo in itertools.combinations_with_replacement(subs,mm):
                        if sum(c[2] for c in combo)==N and not (mm==1 and combo[0][0]=='BA'):
                            jobs.append((comb(v,k)//(q*d),'F',N,t,k,v,n,inf,d,combo))
        elif not isprime(n):   # Z_n composite
            units=[a for a in range(1,n) if gcd(a,n)==1]
            cyc={}
            for a in units:
                S=set(); y=1
                while True:
                    y=y*a%n; S.add(y)
                    if y==1: break
                cyc.setdefault(frozenset(S),a)
            for S,a in cyc.items():
                d=len(S)
                if d<2 or comb(v,k)//(n*d)>CAP: continue
                subs=[]
                for s in range(1,n+1):
                    if n%s: continue
                    for e in range(1,d+1):
                        if d%e==0: subs.append((s,e,s*d//e))
                for mm in range(1,MAXM+1):
                    for combo in itertools.combinations_with_replacement(subs,mm):
                        if sum(c[2] for c in combo)==N:
                            jobs.append((comb(v,k)//(n*d),'Z',N,t,k,v,n,inf,(a,d),combo))
jobs.sort(key=lambda j:j[0]); print(len(jobs),"jobs",len({j[2:6] for j in jobs}),"targets",flush=True)
done=set()
for job in jobs:
    _,typ,N,t,k,v,n,inf,par,combo=job
    if (N,t,k,v) in done: continue
    ext=[n] if inf else []
    if typ=='Z':
        a,d=par
        sig=tuple([(x+1)%n for x in range(n)]+ext)
        def mu(c): return tuple([(c*x)%n for x in range(n)]+ext)
        H=[sig,mu(a)]; Ks=[]
        for s,e,_ in combo:
            g=[]
            if s<n: g.append(tuple([(x+s)%n for x in range(n)]+ext))
            if e>1: g.append(mu(pow(a,d//e,n)))
            Ks.append(g)
    else:
        d=par; p,m=ppow(n)
        if n not in fields: fields[n]=gf(p,m)
        mul,add,w=fields[n]
        def pw(x,e):
            y=1
            for _ in range(e): y=mul(y,x)
            return y
        def tr(b): return tuple([add(x,b) for x in range(n)]+ext)
        def mu(c): return tuple([mul(c,x) for x in range(n)]+ext)
        a=pw(w,(n-1)//d)
        basis=[p**i for i in range(m)]
        H=[tr(b) for b in basis]+[mu(a)]; Ks=[]
        for B,e,_ in combo:
            g=[mu(pw(a,d//e))] if e>1 else []
            if B=='BA': g+= [tr(b) for b in basis]
            elif B.startswith('BF'):
                s=int(B[2:]); z=pw(w,(n-1)//(s-1))   # generator of F_s^*
                g+= [tr(1)] + [tr(pw(z,i)) for i in range(1,s-1)]
            Ks.append(g)
    lam=comb(v-t,k-t)//N; t0=time.time()
    tag=f"LS[{N}]({t},{k},{v}) {typ} X={'Z' if typ=='Z' else 'F'}{n}{'+inf' if inf else ''} par={par} combo={combo}"
    try: st,sol,info=solve(v,t,k,lam,H,Ks,tl=TL)
    except AssertionError as ex: print("SKIP",tag,ex,file=out,flush=True); continue
    if st=="SAT":
        ok=verify(v,t,k,lam,H,Ks,sol); print(tag,"SAT ok=",ok,"H=",H,"Ks=",Ks,"sol=",sol,file=out,flush=True)
        if ok: done.add((N,t,k,v))
    else: print(tag,st,f"[{time.time()-t0:.0f}s]",file=out,flush=True)
