#!/usr/bin/env python3
"""Verify the shift-generated large sets in data/shift_generated_large_sets.txt.

Each entry gives N, t, k, v and lambda, generators of a group G and a shifter sigma, all as image
lists on the points 0,...,v-1 (point i is mapped to the i-th number), and base blocks.  The first
design B is the union of the G-orbits of the base blocks, and the large set consists of the images
of B under sigma^0, sigma^1, ..., sigma^(N-1).

The program generates G from its generators, expands the orbits, and checks that
  (1) B is a t-(v,k,lambda) design, and
  (2) the N images of B under the powers of sigma are pairwise disjoint and together contain
      every k-subset.
Every k-subset and every t-subset is handled explicitly.  Only the Python standard library is used.
Usage:  python3 verify_shift_generated.py [datafile]
"""
import re
import sys
from itertools import combinations
from math import comb


def read(path):
    entries, cur = [], None
    for line in open(path):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        m = re.match(r'LS\[(\d+)\]\((\d+),(\d+),(\d+)\) lambda=(\d+)', line)
        if m:
            N, t, k, v, lam = map(int, m.groups())
            cur = dict(N=N, t=t, k=k, v=v, lam=lam, gens=[], sigma=None, blocks=[], name=line)
            entries.append(cur)
        elif line.startswith('G:'):
            cur['gens'].append(tuple(map(int, line[2:].split())))
        elif line.startswith('sigma:'):
            cur['sigma'] = tuple(map(int, line[6:].split()))
        elif line.startswith('blocks:'):
            cur['blocks'] = [tuple(map(int, b.split())) for b in line[7:].split('|')]
    return entries


def group(gens, v):
    e = tuple(range(v))
    G, frontier = {e}, [e]
    while frontier:
        new = []
        for g in frontier:
            for h in gens:
                gh = tuple(h[g[i]] for i in range(v))
                if gh not in G:
                    G.add(gh)
                    new.append(gh)
        frontier = new
    return G


def check(e):
    N, t, k, v, lam, s = e['N'], e['t'], e['k'], e['v'], e['lam'], e['sigma']
    errors = []
    G = group(e['gens'], v)
    B = set()
    for b in e['blocks']:
        for g in G:
            B.add(frozenset(g[x] for x in b))
    if any(len(K) != k for K in B):
        errors.append('a block has the wrong size')
    count = {}
    for K in B:
        for T in combinations(sorted(K), t):
            count[T] = count.get(T, 0) + 1
    if len(count) != comb(v, t) or any(c != lam for c in count.values()):
        errors.append(f'the first design is not a {t}-({v},{k},{lam}) design')
    seen, image = set(), B
    for i in range(N):
        if seen & image:
            errors.append(f'the image under sigma^{i} meets an earlier design')
        seen |= image
        image = {frozenset(s[x] for x in K) for K in image}
    if len(seen) != comb(v, k):
        errors.append(f'the designs contain {len(seen)} of the {comb(v, k)} {k}-subsets')
    return errors


if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else 'data/shift_generated_large_sets.txt'
    bad = 0
    for e in read(path):
        err = check(e)
        print(e['name'], ': VERIFIED' if not err else ': FAILED', flush=True)
        for m in err:
            print('   ', m)
        bad += bool(err)
    sys.exit(1 if bad else 0)
