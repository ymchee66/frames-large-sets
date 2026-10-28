"""Verify the large sets of data/rank_two_large_sets.json.
For each entry: expand each base design as the union of the orbits of its base blocks under K_i,
take all images under H, and check that the N resulting sets are t-(v,k,lambda) designs that
partition the k-subsets of {0,...,v-1}."""
import json, sys
from itertools import combinations
from math import comb
def close(gens, v):
    e=tuple(range(v)); G={e}; fr=[e]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                c=tuple(h[g[i]] for i in range(v))
                if c not in G: G.add(c); nw.append(c)
        fr=nw
    return G
def check(r):
    N,t,k,v,lam=r['N'],r['t'],r['k'],r['v'],r['lam']
    H=close([tuple(g) for g in r['H_generators']],v)
    designs=set()
    for Kg,base in zip(r['K_generators'],r['base_blocks']):
        K=close([tuple(g) for g in Kg],v)
        B=frozenset(tuple(sorted(g[x] for x in b)) for g in K for b in base)
        for h in H: designs.add(frozenset(tuple(sorted(h[x] for x in b)) for b in B))
    if len(designs)!=N: return False
    allb=[b for D in designs for b in D]
    if len(allb)!=comb(v,k) or len(set(allb))!=comb(v,k): return False
    for D in designs:
        cnt={}
        for b in D:
            for T in combinations(b,t): cnt[T]=cnt.get(T,0)+1
        if len(cnt)!=comb(v,t) or set(cnt.values())!={lam}: return False
    return True
if __name__=="__main__":
    path=sys.argv[1] if len(sys.argv)>1 else "data/rank_two_large_sets.json"
    ok=True
    for r in json.load(open(path)):
        res=check(r); ok&=res
        print(f"LS[{r['N']}]({r['t']},{r['k']},{r['v']}) lambda={r['lam']} orbits={'+'.join(map(str,r['indices']))}: {'OK' if res else 'FAILED'}",flush=True)
    print("all OK" if ok else "FAILURES")
