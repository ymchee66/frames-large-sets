"""Small permutation groups with multiplication table; subgroup classes by bottom-up joins."""
import numpy as np
from itertools import combinations
class Group:
    def __init__(s, gens, v):
        s.v=v; e=np.arange(v,dtype=np.int16)
        idx={e.tobytes():0}; E=[e]; fr=[e]
        G=[np.array(g,dtype=np.int16) for g in gens]
        while fr:
            nw=[]
            for a in fr:
                for g in G:
                    c=g[a]          # apply a then g
                    b=c.tobytes()
                    if b not in idx: idx[b]=len(E); E.append(c); nw.append(c)
            fr=nw
        s.E=np.array(E); s.idx=idx; s.n=len(E)
        n=s.n; M=np.empty((n,n),dtype=np.int32)
        for i in range(n):
            C=s.E[:,s.E[i]]        # row j: apply E[i] then E[j]
            M[i]=[idx[r.tobytes()] for r in C]
        s.M=M                       # M[i,j] = E[i] then E[j]
        s.inv=np.empty(n,dtype=np.int32)
        for i in range(n): s.inv[i]=int(np.where(M[i]==0)[0][0])
        s.gens_idx=[idx[g.tobytes()] for g in G]
    def closure(s, S):
        gens=np.array(sorted(set(S)),dtype=np.int64)
        mask=np.zeros(s.n,dtype=bool); mask[0]=True
        if len(gens)==0: return frozenset([0])
        fr=np.array([0]); 
        while len(fr):
            nw=np.unique(s.M[fr][:,gens].ravel())
            nw=nw[~mask[nw]]
            mask[nw]=True; fr=nw
        return frozenset(np.nonzero(mask)[0].tolist())
    def conj(s, K, h):   # h^-1 K h
        hi=s.inv[h]; return frozenset(int(s.M[s.M[hi,k],h]) for k in K)
    def perm(s,i): return tuple(int(x) for x in s.E[i])
    def gens_of(s,K):
        if hasattr(s,'G_of') and K in s.G_of: return s.G_of[K]
        return s._gens_of(K)
    def _gens_of(s,K):
        """small generating set of subgroup K"""
        gens=[]; cur=frozenset([0])
        for k in sorted(K, key=lambda x:-s.order(x)):
            if k not in cur:
                gens.append(k); cur=s.closure(gens)
                if len(cur)==len(K): break
        return gens
    def order(s,i):
        o=1; a=i
        while a!=0: a=int(s.M[a,i]); o+=1
        return o
    def invariant(s,K):
        ords=sorted(s.order(k) for k in K)
        # orbit lengths on points
        P=s.E[sorted(K)]
        seen=set(); lens=[]
        for x in range(s.v):
            if x in seen: continue
            o=set(P[:,x].tolist()); seen|=o; lens.append(len(o))
        return (len(K),tuple(ords),tuple(sorted(lens)))
    def subgroup_classes(s, max_index=None, max_classes=5000):
        """conjugacy-class representatives of all subgroups (via joins of cyclic subgroups)."""
        cyc={}
        def ppow(o):
            p=next(d for d in range(2,o+1) if o%d==0)
            while o%p==0: o//=p
            return o==1
        for i in range(1,s.n):
            if not ppow(s.order(i)): continue
            C=s.closure([i])
            if C not in cyc: cyc[C]=[i]
        s.G_of=dict(cyc)
        cyc=list(cyc.keys())
        s.cyc_of=[None]*s.n
        for j,C in enumerate(cyc):
            o=len(C)
            for x in C:
                if s.order(x)==o: s.cyc_of[x]=j
        classes={}; seen=set()
        Mi=s.M[s.inv]   # Mi[h,k] = h^-1 * k
        def add(K):
            if K in seen: return False
            classes.setdefault(len(K),[]).append(K)
            Ka=np.array(sorted(K))
            C=s.M[Mi[:,Ka], np.arange(s.n)[:,None]]   # row h: h^-1 k h
            C.sort(axis=1)
            for row in C: seen.add(frozenset(row.tolist()))
            return True
        allsub=[]
        for C in cyc:
            if add(C): allsub.append(C)
        add(frozenset([0])); s.G_of[frozenset([0])]=[]
        cid={C:i for i,C in enumerate(cyc)}
        cgen=np.array([s.G_of[C][0] for C in cyc])
        frontier=list(allsub)
        while frontier:
            nw=[]
            for A in frontier:
                Ka=np.array(sorted(A))
                Cj=s.M[Mi[:,Ka], np.arange(s.n)[:,None]]; Cj.sort(axis=1)
                Nrm=[h for h in range(s.n) if np.array_equal(Cj[h],Ka)]
                visited=np.zeros(len(cyc),dtype=bool)
                Amask=np.zeros(s.n,dtype=bool); Amask[Ka]=True
                for j,C in enumerate(cyc):
                    if visited[j]: continue
                    for h in Nrm:
                        g=int(s.M[s.M[s.inv[h],cgen[j]],h])
                        visited[cid[s.closure([g])] if False else j]=True
                    # mark orbit via generator images
                    imgs={int(s.M[s.M[s.inv[h],cgen[j]],h]) for h in Nrm}
                    for g in imgs:
                        Cg=s.cyc_of[g]
                        if Cg is not None: visited[Cg]=True
                    if Amask[cgen[j]]: continue
                    g=s.G_of[A]+s.G_of[C]
                    B=s.closure(g)
                    if B not in s.G_of: s.G_of[B]=g
                    if add(B): nw.append(B)
            frontier=nw
            if sum(len(L) for L in classes.values())>max_classes: raise RuntimeError("too many")
        out=[R for L in classes.values() for R in L]
        if max_index: out=[R for R in out if s.n//len(R)<=max_index]
        return out
    def is_conj(s,A,B):
        if len(A)!=len(B): return False
        if A==B: return True
        a=min(A-{0}) if len(A)>1 else 0
        for h in range(s.n):
            if s.conj(A,h)==B: return True
        return False
