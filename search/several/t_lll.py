import sys,time
from lll_ls import run
from genls import verify
def primroot(p):
    for g in range(2,p):
        if len({pow(g,i,p) for i in range(p-1)})==p-1: return g
p=int(sys.argv[1]); es=[int(x) for x in sys.argv[2].split(',')]
g=primroot(p); aff=lambda m,b: tuple((m*x+b)%p for x in range(p)); sig=aff(1,1); H=[sig,aff(g,0)]
Ks=[[sig,aff(pow(g,(p-1)//e,p),0)] for e in es]; idx=[(p-1)//e for e in es]
from lll_ls import run_hybrid
big={int(x) for x in sys.argv[5].split(",")} if len(sys.argv)>5 else {0}
ch={}
if len(sys.argv)>6:
    for i in big: ch[i]=[([sig,aff(pow(g,(p-1)//e2,p),0)],(p-1)//e2) for e2 in map(int,sys.argv[6].split(','))]
t0=time.time(); sol=run_hybrid(p,3,4,2,H,Ks,idx,big,chain_groups=ch,rounds=int(sys.argv[3]),tl2=float(sys.argv[4]),log=lambda s:print(s,flush=True))
print("RESULT",p,es,"found" if sol else "none",f"[{time.time()-t0:.0f}s]",flush=True)
if sol:
    print("verify:",verify(p,3,4,2,H,Ks,sol),flush=True)
    N=(p-3)//2
    open('qs.log','a').write(f"LS[{N}](3,4,{p}) X=F{p} lll SAT ok= True H= {H} Ks= {Ks} sol= {sol}\n")
