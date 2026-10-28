import json
def close(gens,v):
    e=tuple(range(v)); G={e}; fr=[e]
    while fr:
        nw=[]
        for g in fr:
            for h in gens:
                c=tuple(h[g[i]] for i in range(v))
                if c not in G: G.add(c); nw.append(c)
        fr=nw
    return G
def inv(h):
    r=[0]*len(h)
    for i,x in enumerate(h): r[x]=i
    return tuple(r)
def comp(a,b): return tuple(b[a[i]] for i in range(len(a)))
def info(r):
    v=r['target'][3]
    H=close([tuple(g) for g in r['Hgens']],v)
    cores=[]
    for Kg in r['Kgens']:
        K=close([tuple(g) for g in Kg],v); c=set(K)
        for h in H:
            hi=inv(h); c&={comp(comp(hi,g),h) for g in K}
        cores.append(c)
    G=cores[0]&cores[1]
    Hl=list(H); ab=None
    # P abelian iff commutators in G
    import itertools, random
    ab=all(comp(comp(comp(inv(a),inv(b)),a),b) in G for a in Hl[:60] for b in Hl[:60])
    return len(H),len(G),len(H)//len(G),ab
if __name__=="__main__":
    R=json.load(open('rank2_realised.json'))
    out={}
    for k,r in R.items():
        out[k]=info(r)
    json.dump(out,open('pinfo.json','w'))
    for k,x in sorted(out.items(),key=lambda z:(json.loads(z[0])[3],json.loads(z[0])[0])):
        r=R[k]; print(json.loads(k),r['H'],r['idx'],x)
