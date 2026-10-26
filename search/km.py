"""Kramer-Mesner matrices for AGL(1,p), ASL(1,p) on GF(p) and PGL(2,p), PSL(2,p) on GF(p)+{inf}, p prime.
Orbits are found through canonical forms obtained from the sharply 2- (resp. 3-) transitive structure."""
from itertools import combinations, permutations
from math import comb
INF = -1
class Grp:
    def __init__(s, name, p):
        s.name, s.p = name, p
        s.inv = [0] + [pow(a, p - 2, p) for a in range(1, p)]
        s.sq = [False] * p
        for a in range(1, p): s.sq[a * a % p] = True
        if name in ("AGL", "ASL"):
            s.v = p; s.pts = list(range(p)); s.base = 2
            s.order = p * (p - 1) // (2 if name == "ASL" else 1)
        else:
            s.v = p + 1; s.pts = [INF] + list(range(p)); s.base = 3
            s.order = p * (p * p - 1) // (2 if name == "PSL" else 1)
    def images(s, S):
        """all images of S under group elements sending an ordered base-tuple of S to the standard base."""
        p, inv = s.p, s.inv; out = []
        if s.base == 2:
            for a, b in permutations(S, 2):
                d = (b - a) % p
                if s.name == "ASL" and not s.sq[d]: continue
                m = inv[d]
                out.append(tuple(sorted((x - a) * m % p for x in S)))
        else:
            for a, b, c in permutations(S, 3):
                # x -> (x-b)(c-a)/((x-a)(c-b)), sending a->inf, b->0, c->1
                if a == INF:   A, B, C, D = 1, -b, 0, (c - b); det = (c - b)
                elif b == INF: A, B, C, D = 0, (c - a), 1, -a; det = -(c - a)
                elif c == INF: A, B, C, D = 1, -b, 1, -a; det = (b - a)
                else: A, B, C, D = (c - a), -b * (c - a), (c - b), -a * (c - b); det = (c - a) * (c - b) * (b - a)
                if s.name == "PSL" and not s.sq[det % p]: continue
                im = []
                for x in S:
                    if x == INF:
                        im.append(INF if C % p == 0 else A * inv[C % p] % p)
                    else:
                        den = (C * x + D) % p
                        im.append(INF if den == 0 else (A * x + B) * inv[den] % p)
                out.append(tuple(sorted(im)))
        return out
    def canon(s, S):
        ims = s.images(S); m = min(ims)
        return m, ims.count(m)          # canonical form, |stabiliser|
    def orbit_reps(s, k):
        base = [0, 1] if s.base == 2 else [INF, 0, 1]
        rest = [x for x in s.pts if x not in base]
        reps = {}
        for R in combinations(rest, k - s.base):
            c, st = s.canon(tuple(base) + R)
            if c not in reps: reps[c] = s.order // st
        assert sum(reps.values()) == comb(s.v, k), (sum(reps.values()), comb(s.v, k))
        return reps                      # canonical rep -> orbit length
    def km(s, t, k):
        rows = s.orbit_reps(t); cols = s.orbit_reps(k)
        cidx = {c: j for j, c in enumerate(sorted(cols))}
        A = []
        for T in sorted(rows):
            r = [0] * len(cidx); rest = [x for x in s.pts if x not in T]
            for R in combinations(rest, k - t):
                c, _ = s.canon(T + R); r[cidx[c]] += 1
            A.append(r)
        return sorted(rows), sorted(cols), [cols[c] for c in sorted(cols)], A

class GrpInf(Grp):
    """AGL(1,p) or ASL(1,p) acting on GF(p)+{inf}, fixing inf (a point stabiliser in PGL(2,p), PSL(2,p))."""
    def __init__(s, name, p):
        Grp.__init__(s, name[:3], p); s.fullname = name
        s.v = p + 1; s.pts = [INF] + list(range(p))
    def canon(s, S):
        if INF in S:
            c, st = Grp.canon(s, tuple(x for x in S if x != INF)); return (INF,) + c, st
        return Grp.canon(s, S)
    def orbit_reps(s, k):
        reps = {}
        for kk, add in ((k, ()), (k - 1, (INF,))):
            rest = list(range(2, s.p))
            for R in combinations(rest, kk - 2):
                c, st = Grp.canon(s, (0, 1) + R)
                if add + c not in reps: reps[add + c] = s.order // st
        assert sum(reps.values()) == comb(s.v, k)
        return reps
