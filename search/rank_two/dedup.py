"""Group phase-1 candidates by an invariant of the frame (orders, point orbits, cycle-type distributions)."""
import sys, json
from collections import Counter
import numpy as np
def elements(gens,v):
    e=np.arange(v,dtype=np.int16); seen={e.tobytes()}; E=[e]; fr=[e]
    G=[np.array(g,dtype=np.int16) for g in gens]
    while fr:
        nw=[]
        for a in fr:
            for g in G:
                c=g[a]; b=c.tobytes()
                if b not in seen: seen.add(b); E.append(c); nw.append(c)
        fr=nw
    return np.array(E)
def ctype(p):
    v=len(p); seen=np.zeros(v,bool); L=[]
    for x in range(v):
        if not seen[x]:
            l=0;y=x
            while not seen[y]: seen[y]=True; y=p[y]; l+=1
            L.append(l)
    return tuple(sorted(L))
CACHE={}
def inv(gens,v):
    key=json.dumps(gens)
    if key not in CACHE:
        E=elements(gens,v)
        orbs=[]; seen=set()
        for x in range(v):
            if x in seen: continue
            o=set(E[:,x].tolist()); seen|=o; orbs.append(len(o))
        CACHE[key]=(len(E),tuple(sorted(orbs)),tuple(sorted(Counter(ctype(p) for p in E).items())))
    return CACHE[key]
classes={}
for l in open(sys.argv[1]):
    r=json.loads(l); v=r['target'][3]
    key=(tuple(r['target']),inv(r['Hgens'],v),tuple(inv(K,v) for K in r['Kgens']))
    classes.setdefault(key,[]).append(r)
out=open(sys.argv[2],'w')
for key,L in classes.items():
    r=L[0]; r['nequiv']=len(L); print(json.dumps(r),file=out)
print(len(classes),"classes from",sum(len(L) for L in classes.values()))
