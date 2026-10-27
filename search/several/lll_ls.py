"""Two-stage heuristic for Theorem 4.1:
 stage 1 (N4): partition H-orbits of k-sets into unions U_i (H-invariant t-(v,k,n_i*lam) designs), CP-SAT with random seeds;
 stage 2: split each U_i into one K_i-invariant design whose H-images tile U_i: 0/1 system A x = b solved by
          Kreher-Radziszowski / Wassermann lattice basis reduction (LLL + BKZ)."""
import random, time
from itertools import combinations
try:   # optional: lattice basis reduction (only used when USE_LLL is True)
    from fpylll import IntegerMatrix, LLL, BKZ
except ImportError:
    IntegerMatrix = LLL = BKZ = None
from ortools.sat.python import cp_model
from fix1 import orbits

def lll01(A, b, tries=4, blocks=(20,30,40,50), seed=0, tl=600):
    """Kreher-Radziszowski / Wassermann: find x in {0,1}^n with A x = b.
    Embed, LLL, keep the kernel sublattice (tail = 0), BKZ 2.0 with pruning strategies, scan, then pruned enumeration."""
    from fpylll import GSO, Enumeration, EnumerationError
    m=len(A); n=len(A[0]); rng=random.Random(seed); t0=time.time()
    Nw=max(4*n,100)
    def decode(vec,perm):
        if vec[n] not in (1,-1) or not all(x in (1,-1) for x in vec[:n]): return None
        s=vec[n]
        for cand in ([(1-s*x)//2 for x in vec[:n]],[(1+s*x)//2 for x in vec[:n]]):
            x=[0]*n
            for r2,j in enumerate(perm): x[j]=cand[r2]
            if all(sum(A[i][j]*x[j] for j in range(n) if A[i][j])==b[i] for i in range(m)): return x
        return None
    for attempt in range(tries):
        perm=list(range(n)); rng.shuffle(perm)
        B=IntegerMatrix(n+1, n+1+m)
        for r,j in enumerate(perm):
            B[r,r]=2
            for i in range(m): B[r,n+1+i]=Nw*A[i][j]
        for r in range(n): B[n,r]=1
        B[n,n]=1
        for i in range(m): B[n,n+1+i]=Nw*b[i]
        LLL.reduction(B)
        ker=[[B[r,c] for c in range(n+1)] for r in range(B.nrows) if not any(B[r,c] for c in range(n+1,B.ncols))]
        for vec in ker:
            x=decode(vec,perm)
            if x: return x
        if not ker: return None
        K=IntegerMatrix.from_matrix(ker); LLL.reduction(K)
        def scanK():
            for r in range(K.nrows):
                x=decode([K[r,c] for c in range(K.ncols)],perm)
                if x: return x
            return None
        for bs in blocks:
            if time.time()-t0>tl: return None
            par=BKZ.Param(block_size=min(bs,K.nrows), strategies=BKZ.DEFAULT_STRATEGY, max_loops=4,
                          flags=BKZ.AUTO_ABORT|BKZ.MAX_LOOPS|BKZ.GH_BND)
            BKZ.reduction(K,par)
            x=scanK()
            if x: return x
        # pruned enumeration for vectors of squared norm n+1 (these are exactly the +-1 solutions when middle coord is +-1)
        try:
            M=GSO.Mat(K); M.update_gso()
            enum=Enumeration(M, nr_solutions=200, max_aux_sols=0)
            from fpylll import Pruning
            sols=enum.enumerate(0, K.nrows, n+1, 0, pruning=None) if K.nrows<=60 else []
        except Exception:
            sols=[]
        for _,coef in sols:
            vec=[sum(int(round(coef[i]))*K[i,c] for i in range(K.nrows)) for c in range(K.ncols)]
            x=decode(vec,perm)
            if x: return x
    return None

def setup(v,t,k,Hgens,Ks,idxs):
    kH,kHr=orbits(v,k,Hgens); Hsize={}
    for S,d in kH.items(): Hsize[d]=Hsize.get(d,0)+1
    tH,tHr=orbits(v,t,Hgens)
    AH=[]
    for T in tHr:
        row={}
        for R in combinations([z for z in range(v) if z not in T],k-t):
            d=kH[tuple(sorted(T+R))]; row[d]=row.get(d,0)+1
        AH.append(row)
    types=[]
    for i,K in enumerate(Ks):
        kK,kKr=orbits(v,k,K); Ksz={}
        for S,j in kK.items(): Ksz[j]=Ksz.get(j,0)+1
        tK,tKr=orbits(v,t,K)
        allowed={}   # H-orbit d -> list of allowed K-orbits j
        for j,S in enumerate(kKr):
            d=kH[S]
            if idxs[i]*Ksz[j]==Hsize[d]: allowed.setdefault(d,[]).append(j)
        AK=[]
        for T in tKr:
            row={}
            for R in combinations([z for z in range(v) if z not in T],k-t):
                j=kK[tuple(sorted(T+R))]; row[j]=row.get(j,0)+1
            AK.append(row)
        types.append(dict(kKr=kKr,allowed=allowed,AK=AK))
    return dict(kHr=kHr,AH=AH,types=types)

def stage1(S,idxs,lam,seed,tl=60):
    m=cp_model.CpModel(); x={}; nH=len(S['kHr'])
    for i,T in enumerate(S['types']):
        for d in T['allowed']: x[i,d]=m.NewBoolVar("")
    for d in range(nH):
        c=[x[i,d] for i in range(len(idxs)) if (i,d) in x]
        if not c: return None
        m.AddExactlyOne(c)
    for i in range(len(idxs)):
        for row in S['AH']: m.Add(sum(a*x[i,d] for d,a in row.items() if (i,d) in x)==idxs[i]*lam)
    rng=random.Random(seed)
    m.Maximize(sum(rng.randint(0,1000)*vv for vv in x.values()))   # random objective for diversity
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=1; s.parameters.random_seed=seed
    st=s.Solve(m)
    if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return "INFEASIBLE" if st==cp_model.INFEASIBLE else None
    return [[d for d in range(nH) if (i,d) in x and s.Value(x[i,d])] for i in range(len(idxs))]

SMALL=10**9
STAGE2_TL=120
USE_LLL=False
def stage2(T,U,lam,seed,tl):
    vars_=[j for d in U for j in T['allowed'][d]]; col={j:c for c,j in enumerate(vars_)}
    rows=[]; rhs=[]
    for d in U:
        r=[0]*len(vars_)
        for j in T['allowed'][d]: r[col[j]]=1
        rows.append(r); rhs.append(1)
    for row in T['AK']:
        r=[0]*len(vars_); nz=False
        for j,a in row.items():
            if j in col: r[col[j]]=a; nz=True
        if not nz:
            if lam: return None
            continue
        rows.append(r); rhs.append(lam)
    n=len(vars_)
    x=None
    if n<=SMALL:   # exact check for small systems (fast rejection of a bad stage-1 partition)
        m=cp_model.CpModel(); X=[m.NewBoolVar("") for _ in range(n)]
        for r,b in zip(rows,rhs): m.Add(sum(a*X[c] for c,a in enumerate(r) if a)==b)
        sv=cp_model.CpSolver(); sv.parameters.max_time_in_seconds=STAGE2_TL; sv.parameters.num_workers=2
        st=sv.Solve(m)
        if st==cp_model.INFEASIBLE: return None
        if st in (cp_model.OPTIMAL,cp_model.FEASIBLE): x=[sv.Value(z) for z in X]
    if x is None and USE_LLL: x=lll01(rows,rhs,seed=seed,tl=tl)
    if x is None: return None
    return [T['kKr'][vars_[c]] for c in range(len(vars_)) if x[c]]

def run(v,t,k,lam,Hgens,Ks,idxs,rounds=20,tl2=300,log=print):
    t0=time.time(); S=setup(v,t,k,Hgens,Ks,idxs)
    log(f"setup {time.time()-t0:.0f}s, H-orbits {len(S['kHr'])}, stage-2 sizes ~ "+str([sum(len(a) for a in T['allowed'].values()) for T in S['types']]))
    for r in range(rounds):
        U=stage1(S,idxs,lam,seed=r)
        if U=="INFEASIBLE": log("stage 1 infeasible (N4)"); return None
        if U is None: log("stage 1 timeout"); continue
        sol=[None]*len(idxs)
        order=sorted(range(len(idxs)),key=lambda i: sum(len(S['types'][i]['allowed'][d]) for d in U[i]))
        for i in order:
            T=S['types'][i]
            t1=time.time(); B=stage2(T,U[i],lam,seed=r,tl=tl2)
            log(f" round {r} type {i} (index {idxs[i]}, {sum(len(T['allowed'][d]) for d in U[i])} vars): {'ok' if B else 'fail'} [{time.time()-t1:.0f}s]")
            if B is None: break
            sol[i]=B
        else:
            return sol
    return None

def stage1_hybrid(S,idxs,lam,big,seed,tl=120,chains=None):
    """CP-SAT: H-level partition for all types; full K-level system for types not in `big`."""
    m=cp_model.CpModel(); x={}; y={}; nH=len(S['kHr'])
    for i,T in enumerate(S['types']):
        for d in T['allowed']: x[i,d]=m.NewBoolVar("")
    for d in range(nH):
        c=[x[i,d] for i in range(len(idxs)) if (i,d) in x]
        if not c: return "INFEASIBLE"
        m.AddExactlyOne(c)
    for i in range(len(idxs)):
        for row in S['AH']: m.Add(sum(a*x[i,d] for d,a in row.items() if (i,d) in x)==idxs[i]*lam)
        if i in big:
            for C in (chains or {}).get(i,[]):   # intermediate level L: L-invariant designs of index [L:K]*lam tiling U_i
                z={}
                for d,js in C['allowed'].items():
                    if (i,d) not in x: continue
                    for j in js: z[j]=m.NewBoolVar("")
                    m.Add(sum(z[j] for j in js)==x[i,d])
                for row in C['AK']: m.Add(sum(a*z[j] for j,a in row.items() if j in z)==C['mult']*lam)
            continue
        T=S['types'][i]
        for d,js in T['allowed'].items():
            for j in js: y[i,j]=m.NewBoolVar("")
            m.Add(sum(y[i,j] for j in js)==x[i,d])
        for row in T['AK']: m.Add(sum(a*y[i,j] for j,a in row.items() if (i,j) in y)==lam)
    rng=random.Random(seed)
    m.Maximize(sum(rng.randint(0,1000)*vv for vv in x.values()))
    s=cp_model.CpSolver(); s.parameters.max_time_in_seconds=tl; s.parameters.num_workers=2; s.parameters.random_seed=seed
    s.parameters.stop_after_first_solution=True
    st=s.Solve(m)
    if st==cp_model.INFEASIBLE: return "INFEASIBLE"
    if st not in (cp_model.OPTIMAL,cp_model.FEASIBLE): return None
    U=[[d for d in range(nH) if (i,d) in x and s.Value(x[i,d])] for i in range(len(idxs))]
    small={i:[S['types'][i]['kKr'][j] for (ii,j) in y if ii==i and s.Value(y[ii,j])] for i in range(len(idxs)) if i not in big}
    return U,small

def run_hybrid(v,t,k,lam,Hgens,Ks,idxs,big,rounds=50,tl2=300,log=print,chain_groups=None):
    t0=time.time(); S=setup(v,t,k,Hgens,Ks,idxs)
    chains={}
    for i,lst in (chain_groups or {}).items():
        chains[i]=[]
        for Lgens,Lidx in lst:   # Lidx=[H:L]; [L:K_i]=idxs[i]/Lidx
            CS=setup(v,t,k,Hgens,[Lgens],[Lidx])['types'][0]
            CS['mult']=idxs[i]//Lidx; chains[i].append(CS)
    log(f"setup {time.time()-t0:.0f}s, H-orbits {len(S['kHr'])}, K-level sizes "+str([sum(len(a) for a in T['allowed'].values()) for T in S['types']])+f", LLL types {sorted(big)}")
    for r in range(rounds):
        t1=time.time(); R=stage1_hybrid(S,idxs,lam,big,seed=r,chains=chains)
        if R=="INFEASIBLE": log("stage 1 infeasible"); return None
        if R is None: log(f" round {r}: stage 1 timeout"); continue
        U,small=R; sol=[small.get(i) for i in range(len(idxs))]
        log(f" round {r}: stage 1 ok [{time.time()-t1:.0f}s]")
        okall=True
        for i in sorted(big):
            t1=time.time(); B=stage2(S['types'][i],U[i],lam,seed=r,tl=tl2)
            log(f"   type {i} (index {idxs[i]}, {sum(len(S['types'][i]['allowed'][d]) for d in U[i])} vars) LLL: {'ok' if B else 'fail'} [{time.time()-t1:.0f}s]")
            if B is None: okall=False; break
            sol[i]=B
        if okall: return sol
    return None
