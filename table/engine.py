"""Existence table for LS[p](t,k,v), p prime: closure of known base results under the known recursive theorems.
Each entry stores a justification string."""
from math import comb
from collections import defaultdict
VMAX = 100; TMAX = 6
PRIMES = [2,3,5,7,11,13,17,19,23,29,31,37,41,43]
def adm(p, t, k, v):
    return 0 <= t <= k <= v and all(comb(v - i, k - i) % p == 0 for i in range(t + 1))
def lam_star(t, k, v):
    l = 1
    while not all((l * comb(v - i, t - i)) % comb(k - i, t - i) == 0 for i in range(t)): l += 1
    return l
def prime_divs(n): return [p for p in PRIMES if n % p == 0]

class Table:
    def __init__(s, p):
        s.p = p; s.E = {}; s.NO = {}; s.queue = []; s.by = defaultdict(set)
    def add(s, t, k, v, why):
        if t < 2 or t > TMAX or k <= t or v > VMAX or v - k < t or k >= v: return
        if (t, k, v) in s.E: return
        assert adm(s.p, t, k, v), ("inadmissible", s.p, t, k, v, why)
        assert (t, k, v) not in s.NO, ("conflict", s.p, t, k, v, why, s.NO.get((t, k, v)))
        s.E[(t, k, v)] = why; s.queue.append((t, k, v)); s.by[(t, k)].add(v)
    def has(s, t, k, v):
        if not (0 < k < v): return False
        if t <= 1: return adm(s.p, t, k, v)          # t = 0 trivial, t = 1 Baranyai and Hartman
        return (t, k, v) in s.E
    def local(s, t, k, v):
        p = s.p
        for i in range(0, t + 1):                     # derived and residual large sets
            for l in range(0, i + 1):
                for j in range(0, l + 1):
                    if i: s.add(t - i, k - j, v - l, f"from ({t},{k},{v}) by derivation")
        s.add(t, v - k, v, f"complement of ({t},{k},{v})")
        if s.has(t, k + 1, v): s.add(t, k + 1, v + 1, f"({t},{k},{v}) and ({t},{k+1},{v})")
        if s.has(t, k - 1, v): s.add(t, k, v + 1, f"({t},{k-1},{v}) and ({t},{k},{v})")
        s.prime_ext(t, k, v)
        if v == 2 * k + 1 and (t % 2 == 0 or p == 2): s.add(t + 1, k + 1, 2 * k + 2, f"Alltop from ({t},{k},{v})")
    def prime_ext(s, t, k, v):
        p = s.p; V = v + 1
        if p * V - p > VMAX + p: return
        for i in range(1, p):                          # Ajoodani-Namini
            for j in range(0, i): s.add(t + 1, p * k + i, p * V + j, f"Ajoodani-Namini from ({t},{k},{v})")
        if t >= 2:
            for i in range(-p + 1, p):
                for j in range(-p, i): s.add(t, p * k + i, p * V + j, f"prime-size extension of ({t},{k},{v})")
    def close(s):
        p = s.p
        for v in range(2, VMAX // p + 3):              # consequences of the large sets of 1-designs
            for k in range(1, v):
                if adm(p, 1, k, v) and k < v: s.prime_ext(1, k, v)
        while True:
            while s.queue: s.local(*s.queue.pop())
            before = len(s.E)
            # Ajoodani-Namini and Khosrovshahi: (t,i,v) for t<i<=k and (t,k,u) give (t,k,u+l(v-t))
            for t in range(2, TMAX + 1):
                for v in range(t + 2, VMAX + 1):
                    k = t + 1
                    while s.has(t, k, v):
                        for u in list(s.by[(t, k)]):
                            w = u + (v - t)
                            while w <= VMAX: s.add(t, k, w, f"({t},{k},{u}) and ({t},i,{v}) for i<={k}"); w += v - t
                        k += 1
            # doubling (Tran van Trung 2019, Corollary 4.5(ii)): (t,i,v) for t<i<=k and k>=2t+1 give (t,k,2v)
            for t in range(2, TMAX + 1):
                for v in range(t + 2, VMAX // 2 + 1):
                    k = t + 1
                    while k < v and s.has(t, k, v):
                        if k >= 2 * t + 1: s.add(t, k, 2 * v, f"doubling of ({t},i,{v}) for i<={k}")
                        k += 1
            # root cases (Khosrovshahi and Tayfeh-Rezaie)
            f = 1
            while p ** f + 2 <= VMAX:
                q = p ** f
                for t in range(2, TMAX + 1):
                    if q + t > VMAX: break
                    kmax = t
                    while s.has(t, kmax + 1, q + t) and kmax + 1 <= (q + t) // 2: kmax += 1
                    top = VMAX if kmax >= (q + t) // 2 else kmax
                    for k in range(t + 1, top + 1):
                        for v in range(k + 1, VMAX + 1):
                            if t <= v % q < k % q: s.add(t, k, v, f"root cases at v={q+t}")
                f += 1
            if len(s.E) == before and not s.queue: break
    def close_no(s):
        ch = True
        while ch:
            ch = False
            for t in range(2, TMAX + 1):
                for v in range(t + 2, VMAX + 1):
                    for k in range(t + 1, v - t + 1):
                        if (t, k, v) in s.NO or not adm(s.p, t, k, v): continue
                        for (a, b, c) in [(t - 1, k - 1, v - 1), (t - 1, k, v - 1), (t, v - k, v)]:
                            if (a, b, c) in s.NO:
                                assert (t, k, v) not in s.E, ("conflict", s.p, t, k, v)
                                s.NO[(t, k, v)] = f"no ({a},{b},{c})"; ch = True; break
