import sys, json, time, os
from cyc2 import *
from run3 import groups_for
res = []
tl = float(sys.argv[1]); vmax = int(sys.argv[2])
jobs = []
for r in (5, 7, 11, 13, 3):
    for t in (2, 3, 4):
        for k in range(t + 1, 8):
            for v in range(2 * k, vmax + 1):
                if any(comb(v - i, k - i) % r for i in range(t + 1)): continue
                if not (v % r < k % r): continue
                jobs.append((comb(v, k), r, t, k, v))
jobs.sort()
print(len(jobs), "candidate parameter sets", flush=True)
for (_, r, t, k, v) in jobs:
    c = v // r; sg = shifter(v, r, c); ls = lam_star(t, k, v); lam = comb(v - t, k - t) // r
    for name, gens in groups_for(v, v % r, r, c) + [("trivial", [])]:
        gens = [g for g in gens if g != tuple(range(v))]
        order = len(group(gens)) if gens else 1
        if comb(v, k) / order > 6000 or comb(v, k) > 400000: continue
        t0 = time.time(); st, sol, nr, nc = cyc_search(v, t, k, r, gens, sg, tl=tl)
        ok = verify(v, t, k, r, gens, sg, sol) if sol else None
        print(f"LS[{r}]({t},{k},{v}) lam={lam}{'*' if lam==ls else ''} G={name} |G|={order} rows={nr} cols={nc} {st} verified={ok} [{time.time()-t0:.0f}s]", flush=True)
        if sol:
            res.append(dict(t=t, k=k, v=v, lam=lam, minimal=(lam == ls), N=r, G=name, order=order, gens=[list(g) for g in gens], sigma=list(sg), reps=sol)); json.dump(res, open("cyc4.json", "w")); break
