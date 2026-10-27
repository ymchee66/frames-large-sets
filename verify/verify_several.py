"""Verify the large sets in data/several_base_large_sets.txt.

For each record the program
  1. closes the generators of each stabiliser K_i to a group and expands the base design B_i as the
     union of the K_i-orbits of its base blocks,
  2. forms the orbit of B_i under H by applying the generators of H to whole designs,
  3. checks that the number of designs is N, that every design is a t-(v,k,lambda) design with
     lambda = C(v-t,k-t)/N, and that the designs are pairwise disjoint and cover all k-subsets.

Usage: python3 verify_several.py [data/several_base_large_sets.txt]
"""
import sys
from itertools import combinations
from math import comb


def read(path):
    recs, cur, mode = [], None, None
    for line in open(path):
        line = line.split("#")[0].strip()
        if not line:
            continue
        w = line.split()
        if w[0] == "LS":
            cur = dict(par=tuple(map(int, w[1:5])), H=[], base=[])
            mode = None
        elif w[0] == "H":
            mode = "H"
        elif w[0] == "K":
            cur["base"].append(([], []))
            mode = "K"
        elif w[0] == "BASE":
            mode = "B"
        elif w[0] == "END":
            recs.append(cur)
            cur = None
        elif mode == "H":
            cur["H"].append(tuple(map(int, w)))
        elif mode == "K":
            cur["base"][-1][0].append(tuple(map(int, w)))
        elif mode == "B":
            cur["base"][-1][1].append(tuple(map(int, w)))
    return recs


def closure(gens, v):
    e = tuple(range(v))
    G, frontier = {e}, [e]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                x = tuple(h[g[i]] for i in range(v))
                if x not in G:
                    G.add(x)
                    new.append(x)
        frontier = new
    return G


def image(design, g):
    return frozenset(tuple(sorted(g[x] for x in b)) for b in design)


def check(rec):
    N, t, k, v = rec["par"]
    lam = comb(v - t, k - t) // N
    designs = set()
    for Kgens, reps in rec["base"]:
        K = closure(Kgens, v)
        B = frozenset(tuple(sorted(g[x] for x in r)) for g in K for r in reps)
        orbit, frontier = {B}, [B]
        while frontier:
            new = []
            for D in frontier:
                for h in rec["H"]:
                    E = image(D, h)
                    if E not in orbit:
                        orbit.add(E)
                        new.append(E)
            frontier = new
        designs |= orbit
    ok = len(designs) == N
    seen = set()
    for D in designs:
        if len(D) != comb(v, k) // N or seen & D:
            ok = False
            break
        seen |= D
        cnt = {}
        for b in D:
            for T in combinations(b, t):
                cnt[T] = cnt.get(T, 0) + 1
        if len(cnt) != comb(v, t) or set(cnt.values()) != {lam}:
            ok = False
            break
    ok = ok and len(seen) == comb(v, k)
    return ok, len(designs), lam


if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "data/several_base_large_sets.txt"
    allok = True
    for rec in read(path):
        ok, n, lam = check(rec)
        N, t, k, v = rec["par"]
        print(f"LS[{N}]({t},{k},{v})  lambda={lam}  designs={n}  base designs={len(rec['base'])}  {'OK' if ok else 'FAIL'}", flush=True)
        allok &= ok
    print("all verified" if allok else "FAILURES")
