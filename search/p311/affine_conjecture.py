"""Exact orbit experiments for Conjecture 3.6; no claim of a general proof."""
from itertools import combinations, permutations
from collections import Counter
import json
import sys
import random
import time
from ortools.sat.python import cp_model


def orbit_data(p):
    inv = [0] + [pow(a, -1, p) for a in range(1, p)]
    def canon(K):
        return min(tuple(sorted((x-a)*inv[(b-a)%p] % p for x in K))
                   for a,b in permutations(K,2))
    def stab(K):
        return sum(tuple(sorted((x-a)*inv[(b-a)%p] % p for x in K)) == K
                   for a,b in permutations(K,2))
    ts = sorted({canon((0,1,x)) for x in range(2,p)})
    ks = sorted({canon((0,1,a,b)) for a,b in combinations(range(2,p),2)})
    ti = {T:i for i,T in enumerate(ts)}
    ss = [stab(K) for K in ks]
    A = [[0]*len(ks) for _ in ts]
    for j,K in enumerate(ks):
        counts = Counter(canon(T) for T in combinations(K,3))
        for T,n in counts.items():
            numerator = n*stab(T)
            assert numerator % ss[j] == 0
            A[ti[T]][j] = numerator//ss[j]
    N=(p-3)//4
    assert len(ts)==(p+1)//6 and len(ks)==N*len(ts)
    assert all(sum(row)==p-3 for row in A)
    shorts=[j for j,s in enumerate(ss) if s==2]
    assert len(shorts)==N
    mean=[j for j,K in enumerate(ks) if sum(K)*inv[4]%p in K]
    starter=mean+[ks.index(canon((0,1,2,3)))]
    assert all(sum(row[j] for j in starter)==4 for row in A)
    allowed=[]
    for K in ks:
        a,b,c,d=K
        colors=set()
        for x,y,z,w in ((a,b,c,d),(a,c,b,d),(a,d,b,c)):
            r=(x-y)*inv[(z-w)%p]%p
            if r in (1,p-1):
                continue
            S=canon((1,p-1,r,(-r)%p))
            colors.add(shorts.index(ks.index(S)))
        allowed.append(sorted(colors))
    return ts,ks,ss,A,shorts,starter,allowed


def solve(p, pairing=False, fixed_mean=False, seconds=30):
    ts,ks,ss,A,shorts,starter,allowed=orbit_data(p)
    N=len(shorts)
    model=cp_model.CpModel()
    colors={j:i for i,j in enumerate(shorts)}
    mean_color=colors[next(j for j in starter if ss[j]==2)]
    x={}
    for j in range(len(ks)):
        options=[colors[j]] if j in colors else (allowed[j] if pairing else range(N))
        if fixed_mean and j in starter:
            options=[i for i in options if i==mean_color]
        for i in options:
            x[j,i]=model.new_bool_var(f'x{j}_{i}')
        model.add(sum(x[j,i] for i in options)==1)
    for row in A:
        for i in range(N):
            model.add(sum(a*x[j,i] for j,a in enumerate(row) if a and (j,i) in x)==4)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=4
    solver.parameters.random_seed=1
    status=solver.solve(model)
    out={'p':p,'pairing':pairing,'fixed_mean':fixed_mean,
         'status':solver.status_name(status),'seconds':solver.wall_time}
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        designs=[[j for j in range(len(ks)) if (j,i) in x and solver.value(x[j,i])]
                 for i in range(N)]
        assert sorted(j for D in designs for j in D)==list(range(len(ks)))
        assert all(sum(row[j] for j in D)==4 for row in A for D in designs)
        out['base_blocks']=[[ks[j] for j in D] for D in designs]
    return out


def peel(p, seconds=150, seed=1):
    """Search with a prescribed mean design; failures imply nothing mathematical."""
    ts,ks,ss,A,shorts,starter,allowed=orbit_data(p)
    rng=random.Random(seed)
    start=time.monotonic()
    designs=[starter]
    shortset=set(shorts)
    failures=0
    deepest=1
    while time.monotonic()-start<seconds:
        used={j for D in designs for j in D}
        remaining=[j for j in range(len(ks)) if j not in used]
        ns=len(shorts)-len(designs)
        if ns==0:
            assert all(sum(row[j] for j in D)==4 for row in A for D in designs)
            return {'p':p,'fixed_mean':True,'status':'VERIFIED_KM',
                    'seconds':time.monotonic()-start,
                    'base_blocks':[[ks[j] for j in D] for D in designs]}
        if ns==1:
            assert all(sum(row[j] for j in remaining)==4 for row in A)
            designs.append(remaining)
            continue
        model=cp_model.CpModel()
        rng.shuffle(remaining)
        x={j:model.new_bool_var(f'x{j}') for j in remaining}
        for row in A:
            model.add(sum(row[j]*x[j] for j in remaining if row[j])==4)
        model.add(sum(x[j] for j in remaining if j in shortset)==1)
        # Random variable order makes successive attempts choose different designs.
        model.add_decision_strategy(list(x.values()),cp_model.CHOOSE_FIRST,
                                    cp_model.SELECT_MAX_VALUE)
        solver=cp_model.CpSolver()
        solver.parameters.max_time_in_seconds=min(2,seconds-(time.monotonic()-start))
        solver.parameters.num_search_workers=1
        solver.parameters.random_seed=rng.randrange(1,1000000)
        solver.parameters.search_branching=cp_model.FIXED_SEARCH
        status=solver.solve(model)
        if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            D=[j for j in remaining if solver.value(x[j])]
            assert all(sum(row[j] for j in D)==4 for row in A)
            designs.append(D)
            if len(designs)>deepest:
                deepest=len(designs)
                if deepest%10==0:
                    print(f'p={p}: {deepest}/{len(shorts)} designs',flush=True)
        else:
            failures+=1
            back=1+min(failures//10,6)
            del designs[max(1,len(designs)-back):]
    return {'p':p,'fixed_mean':True,'status':'TIME_LIMIT',
            'deepest':deepest,'failures':failures,'seconds':time.monotonic()-start}


if __name__=='__main__':
    p=int(sys.argv[1])
    result=solve(p,pairing='--pairing' in sys.argv,fixed_mean='--mean' in sys.argv)
    print(json.dumps(result))
