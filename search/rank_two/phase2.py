"""Phase 2: full filter (N1, N4) and the Theorem 3.1 solver for phase-1 candidates; follows the file."""
import sys, json, time, os, subprocess
from math import comb
from filters import n1_n4
from genls2 import solve, verify
import importlib.util
spec=importlib.util.spec_from_file_location('ddinv','/home/claude/ls/rank2/dedup_inv.py'); dd=importlib.util.module_from_spec(spec); spec.loader.exec_module(dd)
seeninv=set()
cand, outp = sys.argv[1], sys.argv[2]; TL=float(sys.argv[3]) if len(sys.argv)>3 else 300
done=set(); found=set()
if os.path.exists(outp):
    for l in open(outp):
        r=json.loads(l); done.add(r['key'])
        if r['status']=='SAT' and r.get('ok'): found.add(tuple(r['target']))
for key in list(done):
    T_,Hg_,Kg_=json.loads(key)
    seeninv.add((tuple(T_),dd.inv(Hg_,T_[3]),tuple(dd.inv(K,T_[3]) for K in Kg_)))
out=open(outp,'a')
def alive(): return subprocess.run(['pgrep','-f','phase1.py'],capture_output=True).returncode==0
fh=open(cand); buf=''
while True:
    l=fh.readline()
    if not l or not l.endswith('\n'):
        buf+=l
        if not alive() and not l: break
        time.sleep(20); continue
    l=buf+l; buf=''
    r=json.loads(l)
    key=json.dumps([r['target'],r['Hgens'],r['Kgens']])
    if key in done or tuple(r['target']) in found: continue
    done.add(key)
    v_=r['target'][3]
    ik=(tuple(r['target']),dd.inv(r['Hgens'],v_),tuple(dd.inv(K,v_) for K in r['Kgens']))
    if ik in seeninv: continue
    seeninv.add(ik)
    N,t,k,v=r['target']; lam=comb(v-t,k-t)//N
    t0=time.time()
    fl=n1_n4(v,t,k,lam,r['Hgens'],r['Kgens'],r['idx'],tl=60)
    rec=dict(key=key,target=r['target'],H=r['H'],idx=r['idx'],Kord=r['Kord'])
    if fl not in ('pass','N4?'): rec['status']='FILTER '+fl
    else:
        import multiprocessing as mp
        q=mp.Queue()
        def job(q,args):
            try: q.put(solve(*args,tl=TL,workers=1))
            except Exception as e: q.put(("ERR "+str(e),None,{}))
        pr=mp.Process(target=job,args=(q,(v,t,k,lam,r['Hgens'],r['Kgens'],r['idx'])))
        pr.start()
        try: st,sol,info=q.get(timeout=4*TL+120)
        except Exception: st,sol,info="HUNG",None,{}
        pr.join(5)
        if pr.is_alive(): pr.kill()
        rec['status']=st; rec['info']=info
        if st=='SAT':
            rec['ok']=verify(v,t,k,lam,r['Hgens'],r['Kgens'],sol,N); rec['sol']=sol; rec['Hgens']=r['Hgens']; rec['Kgens']=r['Kgens']
            if rec['ok']: found.add(tuple(r['target']))
    rec['time']=round(time.time()-t0)
    print(json.dumps(rec),file=out,flush=True)
