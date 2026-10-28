import json, csv
from math import comb
from collections import defaultdict
T=json.load(open('targets_all20.json'))
st={(int(r['N']),int(r['t']),int(r['k']),int(r['v'])):r['status'] for r in csv.DictReader(open('/home/claude/ls/tab/large_sets_table.csv'))}
res=defaultdict(list)
for g in ('resall.jsonl','resA.jsonl','resB.jsonl','resC.jsonl','resD.jsonl'):
    seen=set()
    for l in open(g):
        r=json.loads(l)
        k=(r['key'],r['status'])
        if k in seen: continue
        seen.add(k); res[tuple(r['target'])].append(r)
cand=defaultdict(int)
for l in open('candall.jsonl'): cand[tuple(json.loads(l)['target'])]+=1
rows=[]
for x in T:
    t4=tuple(x[1:]); L=res[t4]
    sat=[r for r in L if r['status']=='SAT' and r.get('ok')]
    unk=[r for r in L if r['status'] in ('UNKNOWN','HUNG')]
    if sat: o='realised'; ex=sat[0]
    elif cand[t4]==0: o='no frame passes cheap tests'; ex=None
    elif unk: o='undecided'; ex=None
    else: o='all frames excluded'; ex=None
    N,t,k,v=t4
    rows.append(dict(N=N,t=t,k=k,v=v,lam=comb(v-t,k-t)//N,status=st[t4],outcome=o,H=ex['H'] if ex else '',orbits='+'.join(map(str,ex['idx'])) if ex else '',ncand=cand[t4]))
from collections import Counter
print(Counter((r['status'],r['outcome']) for r in rows))
json.dump(rows,open('summary20.json','w'))
