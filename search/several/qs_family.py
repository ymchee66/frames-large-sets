"""LS[(p-3)/2](3,4,p), H=AGL(1,p): enumerate index patterns, prefilter by N1/N2/N4, then solve survivors."""
import sys,time,ast
from genls import solve, verify
from filters import n1_n4
from necessary import count_sub
def primroot(p):
    for g in range(2,p):
        if len({pow(g,i,p) for i in range(p-1)})==p-1: return g
def parts(n,vals,maxparts):
    if n==0: yield []; return
    if maxparts==0: return
    for i,a in enumerate(vals):
        if a<=n:
            for r in parts(n-a,vals[i:],maxparts-1): yield [a]+r
def n2ok(p,es,idx,lam=2,t=3,k=4):
    for q in [q for q in range(2,p) if (p-1)%q==0 and all(q%d for d in range(2,q))]:
        phi=sum(i for e,i in zip(es,idx) if e%q==0); qs=(p-1)//q
        for j in range(t//q+1):
            i=t-j*q
            if i>1 or j>qs: continue
            c=count_sub([1]*(1-i)+[q]*(qs-j),k-t)
            if (c-phi*lam)%q or c<phi*(lam%q): return False
    return True
p=int(sys.argv[1]); TL=float(sys.argv[2]); MAXP=int(sys.argv[3])
N=(p-3)//2; g=primroot(p); aff=lambda m,b: tuple((m*x+b)%p for x in range(p)); sig=aff(1,1); H=[sig,aff(g,0)]
vals=sorted({(p-1)//e for e in range(1,p) if (p-1)%e==0},reverse=True)
surv=[]
for pat in parts(N,vals,MAXP):
    es=[(p-1)//i for i in pat]
    if not n2ok(p,es,pat): print(p,pat,"FILTER N2",flush=True); continue
    Ks=[[sig]+([aff(pow(g,(p-1)//e,p),0)] if e>1 else []) for e in es]
    f=n1_n4(p,3,4,2,H,Ks,pat,tl=60)
    print(p,pat,"FILTER",f,flush=True)
    if f in ("pass","N4?"): surv.append((pat,es,Ks))
print(p,"survivors:",[s[0] for s in surv],flush=True)
for pat,es,Ks in surv:
    t0=time.time(); st,sol,info=solve(p,3,4,2,H,Ks,tl=TL)
    print(p,pat,"e=",es,st,f"[{time.time()-t0:.0f}s]",flush=True)
    if st=="SAT":
        open('qs.log','a').write(f"LS[{N}](3,4,{p}) X=F{p} qs SAT ok= True H= {H} Ks= {Ks} sol= {sol}\n"); break
