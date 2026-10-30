"""Rank-one frames with P = PSL(3,q) acting on the N = q^2+q+1 designs as on the points (or lines) of PG(2,q).
X = PG(2,q) + F (f extra points), q prime."""
import sys, json, time, csv
from math import comb
from itertools import product
from filters import n1_n4
from genls2 import solve, verify
def pg2(q):
    pts=[]
    for x in product(range(q),repeat=3):
        if any(x):
            i=next(j for j in range(3) if x[j]); inv=pow(x[i],q-2,q)
            y=tuple(c*inv%q for c in x)
            if y not in pts: pts.append(y)
    return pts
def norm(x,q):
    i=next(j for j in range(3) if x[j]); inv=pow(x[i],q-2,q); return tuple(c*inv%q for c in x)
def groups(q,f):
    pts=pg2(q); ix={p:i for i,p in enumerate(pts)}; n=len(pts); v=n+f
    w=next(g for g in range(1,q) if q==2 or len({pow(g,e,q) for e in range(q-1)})==q-1)
    def mat(M): return tuple([ix[norm(tuple(sum(M[r][c]*x[c] for c in range(3))%q for r in range(3)),q)] for x in pts]+list(range(n,v)))
    def E(i,j):
        M=[[int(r==c) for c in range(3)] for r in range(3)]; M[i][j]=1; return M
    def D(i):
        M=[[int(r==c) for c in range(3)] for r in range(3)]; M[i][i]=w; return M
    Hg=[mat(E(i,j)) for i in range(3) for j in range(3) if i!=j]+[mat(D(0))]
    # stabiliser of point <e1>: matrices with first column a*e1  (entries M[1][0]=M[2][0]=0)
    Kp=[mat(E(0,1)),mat(E(0,2)),mat(E(1,2)),mat(E(2,1)),mat(D(0)),mat(D(1)),mat(D(2))]
    # stabiliser of line x1=0 (image of <e1> under duality): matrices with first row a*e1
    Kl=[mat(E(1,0)),mat(E(2,0)),mat(E(1,2)),mat(E(2,1)),mat(D(0)),mat(D(1)),mat(D(2))]
    Fsym=[]
    if f>=2:
        Fsym=[tuple(list(range(n))+[n+((i+1)%f) for i in range(f)])]
        if f>=3: Fsym.append(tuple(list(range(n))+[n+1,n]+list(range(n+2,v))))
    return n,v,Hg,{'point':Kp,'line':Kl},Fsym
if __name__=="__main__":
    targets=json.load(open(sys.argv[1])); out=open(sys.argv[2],'a'); TL=float(sys.argv[3]) if len(sys.argv)>3 else 300
    for (N,t,k,v) in targets:
        q={7:2,13:3,31:5}[N]; f=v-N
        n,v2,Hg,Ks,Fsym=groups(q,f); assert v2==v and n==N
        lam=comb(v-t,k-t)//N
        for fp in (['F1','FSym'] if f>=2 else ['F1']):
            for kind in ('point','line'):
                H=Hg+(Fsym if fp=='FSym' else []); K=Ks[kind]+(Fsym if fp=='FSym' else [])
                t0=time.time()
                fl=n1_n4(v,t,k,lam,H,[K],[N],tl=120)
                rec=dict(target=[N,t,k,v],H=f"PSL(3,{q})+{f} {fp}",K=kind,lam=lam)
                if fl not in ('pass','N4?'): rec['status']='FILTER '+fl
                else:
                    st,sol,info=solve(v,t,k,lam,H,[K],[N],tl=TL,workers=1)
                    rec['status']=st; rec['info']=info
                    if st=='SAT':
                        rec['ok']=verify(v,t,k,lam,H,[K],sol,N); rec.update(Hgens=H,Kgens=[K],idx=[N],sol=sol)
                rec['time']=round(time.time()-t0)
                print(json.dumps(rec),file=out,flush=True)
                if rec['status']=='SAT': break
            else: continue
            break
