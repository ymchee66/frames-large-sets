import sys, json, time, os
from cyc2 import *
res = json.load(open("cyc3.json")) if os.path.exists("cyc3.json") else []
def groups_for(v, t, N, c):
    """groups centralising sigma: rotate the c cycles, permute the t fixed points."""
    rot = tuple([(x + N) % (c * N) for x in range(c * N)] + list(range(c * N, v)))
    f = list(range(c * N, v))
    fc = cyc(v, [tuple(f)]); ft = cyc(v, [tuple(f[:2])])
    return [("rot x Sym(fixed)", [rot, fc, ft]), ("rot x cyc(fixed)", [rot, fc]), ("rot", [rot]), ("cyc(fixed)", [fc])]
if __name__ == '__main__':
    tl = float(sys.argv[1])
    for (t, k, v) in [tuple(map(int, a.split(","))) for a in sys.argv[2:]]:
        ls = lam_star(t, k, v); N = (v - t) // ls; sg = shifter(v, N, ls)
        for name, gens in groups_for(v, t, N, ls):
            t0 = time.time(); st, sol, nr, nc = cyc_search(v, t, k, N, gens, sg, tl=tl)
            ok = verify(v, t, k, N, gens, sg, sol) if sol else None
            print(f"({t},{k},{v};{ls}) N={N} G={name} |G|={len(group(gens))} rows={nr} cols={nc} {st} verified={ok} [{time.time()-t0:.0f}s]", flush=True)
            if sol:
                res.append(dict(t=t, k=k, v=v, lam=ls, N=N, G=name, gens=[list(g) for g in gens], sigma=list(sg), reps=sol)); json.dump(res, open("cyc3.json", "w")); break
