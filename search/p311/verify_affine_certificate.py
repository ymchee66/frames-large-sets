"""Verify affine large-set certificates using only Python's standard library.

This uses direct fourth-point completions, not Kramer--Mesner entries or a solver.
Usage: python3 research/verify_affine_certificate.py research/mean_completion_311.json
"""
import json
import sys
from collections import Counter
from itertools import permutations
from math import comb, isqrt


def verify(path):
    with open(path) as f:
        data=json.load(f)
    p=data['p']
    assert p>=11 and p%12==11
    assert all(p%d for d in range(2,isqrt(p)+1)), 'p must be prime'
    N=(p-3)//4
    designs=data['base_blocks']
    assert len(designs)==N
    inv=[0]+[pow(a,-1,p) for a in range(1,p)]
    def key(B):
        return min(tuple(sorted((v-a)*inv[(b-a)%p]%p for v in B))
                   for a,b in permutations(B,2))
    colors={}
    for color,blocks in enumerate(designs):
        assert len(blocks)==(p+1)//6
        short_count=0
        for B in blocks:
            assert len(B)==4 and len(set(B))==4
            assert all(type(x) is int and 0<=x<p for x in B)
            orbit=key(B)
            assert orbit not in colors, 'Repeated affine orbit'
            colors[orbit]=color
            c=sum(B)*inv[4]%p
            short_count+=({(2*c-x)%p for x in B}==set(B))
        assert short_count==1
    # Every triple is an affine image of one tested here. Every four-set
    # has an affine image containing 0 and 1, so coverage is checked too.
    for t in range(2,p):
        count=Counter(colors[key((0,1,t,x))]
                      for x in range(p) if x not in (0,1,t))
        assert len(count)==N and set(count.values())=={4}, (t,count)
    return {'p':p,'designs':N,'affine_orbits':len(colors),
            'normalized_triples_checked':p-2,
            'fourth_point_completions_checked':(p-2)*(p-3),
            'blocks_per_design':comb(p,3),
            'total_blocks':comb(p,4),'verified':True}


if __name__=='__main__':
    print(json.dumps(verify(sys.argv[1]),indent=2))
