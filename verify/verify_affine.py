#!/usr/bin/env python3
"""Verify the large sets in data/affine_large_sets.txt.

Each entry is a large set LS[N](t,k,p) on the points 0,...,p-1 of GF(p), p prime.  Each design is
the union of the orbits of its base blocks under
    AGL(1,p) = {x -> a x + b : a != 0}                 (group AGL), or
    ASL(1,p) = {x -> a x + b : a a nonzero square}      (group ASL, p = 3 mod 4).

The program checks, for every entry, that
  (1) the base blocks of all designs lie in pairwise distinct orbits,
  (2) the orbits together contain every k-subset exactly once, and
  (3) every design contains every t-subset in exactly lambda blocks.

Method.  Both groups act regularly on the set of "admissible" ordered pairs (x,y), x != y
(all pairs for AGL, pairs with y-x a nonzero square for ASL).  For a k-subset K and an admissible
pair (x,y) in K, the unique group element with x -> 0, y -> 1 maps K to a set N(K;x,y).  The
minimum of N(K;x,y) over admissible pairs is a complete orbit invariant, and the number of
admissible pairs attaining the minimum is the order of the stabiliser of K.  Every t-subset with
t >= 2 is the image of a t-subset containing 0 and 1, and the designs are invariant under the
group, so it suffices to check (3) for the t-subsets containing 0 and 1.

Only the Python standard library is used.  Usage:  python3 verify_affine.py [datafile]
"""
import sys
from itertools import combinations
from math import comb


def read(path):
    entries, cur = [], None
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('LS '):
            _, grp, p, t, k, lam, n = line.split()
            cur = dict(group=grp, p=int(p), t=int(t), k=int(k), lam=int(lam), N=int(n), designs=[])
            entries.append(cur)
        else:
            cur['designs'].append([tuple(map(int, b.split())) for b in line.split('|')])
    return entries


def check(e):
    p, t, k, lam, N, grp = e['p'], e['t'], e['k'], e['lam'], e['N'], e['group']
    inv = [0] + [pow(a, p - 2, p) for a in range(1, p)]
    squares = {a * a % p for a in range(1, p)}
    admissible = (lambda d: True) if grp == 'AGL' else (lambda d: d in squares)
    order = p * (p - 1) if grp == 'AGL' else p * (p - 1) // 2

    def canon(K):
        best, stab = None, 0
        for x in K:
            for y in K:
                if x != y and admissible((y - x) % p):
                    m = inv[(y - x) % p]
                    r = tuple(sorted((z - x) * m % p for z in K))
                    if best is None or r < best:
                        best, stab = r, 1
                    elif r == best:
                        stab += 1
        return best, stab

    errors = []
    if len(e['designs']) != N:
        errors.append(f'{len(e["designs"])} designs listed, {N} expected')
    design_of, blocks = {}, 0
    for i, D in enumerate(e['designs']):
        for B in D:
            if len(set(B)) != k or not all(0 <= x < p for x in B):
                errors.append(f'design {i}: bad base block {B}')
                continue
            c, s = canon(B)
            if c in design_of:
                errors.append(f'design {i}: base block {B} repeats an earlier orbit')
            design_of[c] = i
            blocks += order // s
    if blocks != comb(p, k):
        errors.append(f'the orbits contain {blocks} blocks, C({p},{k}) = {comb(p, k)}')
    # design property, for every t-subset containing 0 and 1
    for R in combinations(range(2, p), t - 2):
        T = (0, 1) + R
        count = [0] * N
        rest = [x for x in range(p) if x not in T]
        for S in combinations(rest, k - t):
            c, _ = canon(T + S)
            count[design_of[c]] += 1
        if any(c != lam for c in count):
            errors.append(f'the {t}-subset {T} is in {count} blocks of the designs')
            break
    return errors


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'data/affine_large_sets.txt'
    bad = 0
    for e in read(path):
        err = check(e)
        name = f'LS[{e["N"]}]({e["t"]},{e["k"]},{e["p"]}) under {e["group"]}(1,{e["p"]}), lambda = {e["lam"]}'
        print(name, ': VERIFIED' if not err else ': FAILED', flush=True)
        for m in err:
            print('   ', m)
        bad += bool(err)
    sys.exit(1 if bad else 0)
