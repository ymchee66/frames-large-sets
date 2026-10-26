"""Core utilities: permutation groups on {0..v-1}, orbits on k-subsets, design checks."""
from itertools import combinations
from math import comb, gcd
from collections import Counter

def cyc(v, cycles):
    p = list(range(v))
    for c in cycles:
        for i, x in enumerate(c):
            p[x] = c[(i + 1) % len(c)]
    return tuple(p)

def mul(p, q):            # apply p then q
    return tuple(q[x] for x in p)

def inv(p):
    r = [0] * len(p)
    for i, x in enumerate(p): r[x] = i
    return tuple(r)

def order(p):
    e = tuple(range(len(p))); q = p; n = 1
    while q != e: q = mul(q, p); n += 1
    return n

def group(gens):
    v = len(gens[0]); e = tuple(range(v))
    seen = {e}; frontier = [e]
    while frontier:
        nxt = []
        for g in frontier:
            for h in gens:
                x = mul(g, h)
                if x not in seen: seen.add(x); nxt.append(x)
        frontier = nxt
    return sorted(seen)

def img(S, p):
    return tuple(sorted(p[x] for x in S))

def orbit(S, gens):
    S = tuple(sorted(S)); seen = {S}; fr = [S]
    while fr:
        nx = []
        for T in fr:
            for g in gens:
                U = img(T, g)
                if U not in seen: seen.add(U); nx.append(U)
        fr = nx
    return seen

def orbits(v, k, gens):
    seen = set(); out = []
    for S in combinations(range(v), k):
        if S in seen: continue
        o = orbit(S, gens); seen |= o; out.append(o)
    return out

def tcount(blocks, v, t):
    c = Counter()
    for B in blocks:
        for T in combinations(B, t): c[T] += 1
    vals = set(c.values())
    if len(c) < comb(v, t): vals.add(0)
    return vals

def is_design(blocks, v, t):
    vals = tcount(blocks, v, t)
    return (len(vals) == 1), (next(iter(vals)) if len(vals) == 1 else sorted(vals))

def lam_star(t, k, v):
    l = 1
    while True:
        if all((l * comb(v - i, t - i)) % comb(k - i, t - i) == 0 for i in range(t)): return l
        l += 1

def cycle_type(p):
    seen = [False] * len(p); c = Counter()
    for i in range(len(p)):
        if not seen[i]:
            n = 0; j = i
            while not seen[j]: seen[j] = True; j = p[j]; n += 1
            c[n] += 1
    return tuple(sorted(c.items()))
