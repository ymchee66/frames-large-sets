"""m=1, centraliser family: X = Z_N x Z_c u F, |F| = v mod N, sigma rotates all N-cycles (f fixed points).
H=<A,E>, K=G=<A0,E> with [A:A0]=N; A=<sigma>,A0=1 (diag) or A=Z_N^c (independent rotations), A0 = zero-sum subgroup;
E from: R (rotate the cycles), RF (R times a cycle on F), mu_a (multiplier of order r | N-1), Fc (cycle on F), SF (Sym(F))."""
import sys, csv, time, os
from math import comb, gcd
from genls2 import solve, verify
from filters import n1_n4
CAP=int(sys.argv[1]); TL=float(sys.argv[2]); out=open(sys.argv[3],'a'); TCAP=float(sys.argv[4])
SH,NSH=(int(sys.argv[5]),int(sys.argv[6])) if len(sys.argv)>6 else (0,1)
import genls2
WORKERS=1 if NSH>1 else 2
def gens(N,c,f):
    v=N*c+f; P=lambda f_: tuple(f_(z) for z in range(v))
    def mapxy(fn):
        def g(z):
            if z>=N*c: return z
            y,x=divmod(z,N); x2,y2=fn(x,y); return y2*N+x2
        return g
    sig=P(mapxy(lambda x,y:((x+1)%N,y)))
    tau=[P(mapxy(lambda x,y,yy=yy:((x+1)%N if y==yy else x,y))) for yy in range(c)]
    inv=lambda g: tuple(sorted(range(v),key=lambda i:g[i]))
    comp=lambda a,b: tuple(b[a[i]] for i in range(v))
    zs=[comp(tau[0],inv(tau[y])) for y in range(1,c)]
    R=P(mapxy(lambda x,y:(x,(y+1)%c))) if c>1 else None
    Fc=P(lambda z: z if z<N*c else N*c+((z-N*c+1)%f)) if f>1 else None
    Ft=P(lambda z: {N*c:N*c+1,N*c+1:N*c}.get(z,z)) if f>2 else None
    RF=comp(R,Fc) if (R and Fc) else None
    g0=next((g for g in range(2,N) if len({pow(g,i,N) for i in range(N-1)})==N-1),1)
    mus={}
    for r in range(2,N):
        if (N-1)%r==0:
            a=pow(g0,(N-1)//r,N); mus[r]=P(mapxy(lambda x,y,a=a:(a*x%N,y)))
    return sig,tau,zs,R,Fc,Ft,RF,mus
def frames(N,c,f):
    sig,tau,zs,R,Fc,Ft,RF,mus=gens(N,c,f)
    As=[("diag",[sig],[])]
    if c>1: As.append(("indep",tau,zs))
    Es=[("",[])]
    if R: Es.append(("R",[R]))
    if RF: Es.append(("RF",[RF]))
    if Fc: Es.append(("Fc",[Fc]))
    if R and Fc: Es.append(("R,Fc",[R,Fc]))
    if Ft: Es.append(("SF",[Fc,Ft]))
    if R and Ft: Es.append(("R,SF",[R,Fc,Ft]))
    base=list(Es)
    for r,mu in mus.items():
        for nm,E in base: Es.append((nm+(",mu%d"%r),E+[mu]))
    res=[]
    for an,A,A0 in As:
        for en,E in Es:
            res.append((f"Z{N}x{c}+{f} A={an} E={{{en}}}",A+E,A0+E))
    return res
rows=[r for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','table','existence_table.csv'))) if r['status'] in ('open','exists')]
ST={(int(r['N']),int(r['t']),int(r['k']),int(r['v'])):r['status'] for r in rows}
T=sorted({(comb(int(r['v']),int(r['k'])),int(r['N']),int(r['t']),int(r['k']),int(r['v'])) for r in rows
          if int(r['N'])>=3 and int(r['v'])//int(r['N'])>=2 and int(r['v'])%int(r['N'])<int(r['k'])%int(r['N']) and comb(int(r['v']),int(r['k']))<=CAP})
T=T[SH::NSH]
print(len(T),"targets",flush=True)
prev=set(l.split(' ::')[0] for l in open(sys.argv[3])) if os.path.exists(sys.argv[3]) else set()
for _,N,t,k,v in T:
    lam=comb(v-t,k-t)//N; f=v%N; c=(v-f)//N; t0=time.time(); found=False
    for desc,H,K in sorted([z for z in frames(N,c,f) if 'mu' in z[0]],key=lambda z:-len(z[2])):   # more generators first (larger K)
        if time.time()-t0>TCAP: break
        tag=f"LS[{N}]({t},{k},{v}) [{ST[(N,t,k,v)]}] {desc}"
        if tag in prev: continue
        try:
            fl=n1_n4(v,t,k,lam,H,[K],[N],tl=20)
            if fl not in ("pass","N4?"): print(tag,":: FILTER",fl,file=out,flush=True); continue
            st,sol,info=solve(v,t,k,lam,H,[K],[N],tl=TL,workers=WORKERS)
        except ValueError as ex:
            print(tag,":: ERR",ex,file=out,flush=True); continue
        if st=="SAT":
            ok=verify(v,t,k,lam,H,[K],sol,N)
            print(tag,":: SAT ok=",ok,"H=",H,"Ks=",[K],"sol=",sol,file=out,flush=True)
            if ok: found=True; break
        else: print(tag,":: ",st,info,file=out,flush=True)
    print(f"LS[{N}]({t},{k},{v}) :: DONE found={found} [{time.time()-t0:.0f}s]",file=out,flush=True)
