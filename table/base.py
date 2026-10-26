"""Base data: general theorems and directly constructed large sets, each with its source."""
from math import comb
from engine import *
from itertools import product
def is_pp(q):
    for p in range(2, q + 1):
        if q % p == 0:
            while q % p == 0: q //= p
            return q == 1
# (N, t, k, v, source)
DIRECT = [
 (7,2,4,9,"KMS91"),(14,2,4,10,"KMS91"),(14,2,5,10,"KMS91"),(42,2,5,11,"KMS91"),(15,2,4,12,"KMS91"),(42,2,6,12,"KMS91"),
 (42,3,6,12,"KMS91"),(6,2,5,12,"KMS91"),(6,3,5,12,"KMS91"),
 (55,2,4,13,"Chouinard83"),(3,4,5,13,"KMO93"),(5,3,4,13,"KMO93"),(11,2,4,14,"KMO93"),
 (5,3,4,23,"CCFK90"),(3,4,6,14,"CCFK90"),(3,3,5,13,"CCFK90"),
 (11,2,5,13,"CM98"),(11,2,6,13,"CM98"),(7,3,4,17,"CM98"),
 (2,4,5,12,"Denniston83"),(7,3,5,11,"Eslami03"),
 (5,2,3,17,"LMW01"),(19,2,4,21,"LMW01"),(3,3,7,21,"LMW01"),(3,5,7,24,"LMW01"),(23,2,4,25,"LMW01"),(29,3,5,32,"LMW01"),
 (7,4,5,32,"LMW01"),(3,3,11,30,"LMW01"),(5,3,11,30,"LMW01"),(13,3,11,30,"LMW01"),(23,3,11,30,"LMW01"),
 (39,4,10,33,"LMW01"),(29,4,10,33,"LMW01"),(3,4,6,13,"LMW01"),(7,2,6,17,"LMW01"),(3,2,9,29,"LMW01"),
 (17,2,15,31,"LMW01"),(1885,3,9,32,"LMW01"),(5,3,4,33,"LMW01"),
 (2,3,7,54,"LMW01"),(5,3,7,54,"LMW01"),(7,3,7,54,"LMW01"),(17,3,7,54,"LMW01"),
 (3,5,6,23,"LOTW07"),(3,5,8,24,"LOTW07"),(3,5,7,33,"LOTW07"),(3,6,10,33,"LOTW07"),
 (5,2,8,17,"LOTW07"),(5,3,9,23,"LOTW07"),(5,3,6,28,"LOTW07"),(5,3,9,33,"LOTW07"),(7,3,5,17,"LOTW07"),(7,3,6,17,"LOTW07"),
 (7,3,4,24,"LOTW07"),(7,3,6,24,"LOTW07"),(7,3,13,33,"LOTW07"),(13,2,4,15,"LOTW07"),(13,2,7,15,"LOTW07"),(13,3,12,33,"LOTW07"),
 (13,6,10,33,"LOTW07"),(17,2,4,19,"LOTW07"),(17,2,6,19,"LOTW07"),(17,2,7,19,"LOTW07"),(29,2,7,31,"LOTW07"),(29,3,6,32,"LOTW07"),
 (221,2,8,19,"LOTW07"),(21,3,12,32,"LOTW07"),(29,2,4,31,"LMW01"),
 (3,3,27,84,"LOTW07"),(3,3,28,84,"LOTW07"),(3,3,29,84,"LOTW07"),(3,3,30,84,"LOTW07"),
 (13,2,4,16,"KT04"),(7,3,5,18,"KT04"),(7,3,6,18,"KT04"),(5,3,4,18,"Teirlinck89"),(5,2,4,17,"Teirlinck89"),
 (5,2,5,27,"LMW01"),(5,3,5,28,"LMW01"),
 (5655,3,9,32,"CuMa99"),(70035,3,13,32,"CuMa99"),
 (3,3,8,32,"EN18"),(5,3,8,32,"EN18"),(13,3,8,32,"EN18"),(9,2,4,29,"EN10"),(9,2,5,29,"EN14"),
]
# Table II of LOTW07: q -> (k; primes) for LS[N](3,k,q+1)
T2 = {11:{6:[2]},19:{6:[2],7:[2],8:[2]},
 23:{5:[7],7:[3,19],8:[3],9:[2],10:[2,17,19],12:[2]},25:{12:[2]},
 27:{5:[2,5],6:[2],8:[5],11:[5,19,23],12:[5]},
 29:{10:[3],11:[3,5,13,23],12:[3,5],13:[3,5]},
 31:{8:[3],9:[3,5,13,29],10:[3,13,23,29],11:[3,5,13,23,29],12:[3,7,13,23,29],13:[3,5,7,23,29],14:[3,5,19,23,29],15:[3,19,23,29],16:[3]},
 37:{7:[2],8:[2],12:[2,3,7,17,29,31],13:[2,3,7,17,29,31],14:[2,3,5,17],15:[2,3,5,17,29,31],16:[2,3,5,17,23,29,31],17:[2,3,5,11,23,29,31],18:[2,3,5,7,11,23,29,31]},
 41:{13:[2,17,19,31],15:[2],16:[2],17:[2],19:[2,3,5,13,29,31],20:[2]},
 43:{6:[2],10:[13,19],12:[13,17,19],13:[2,17,19,41],16:[2],17:[2,19,29,31,41],18:[2,3,19,29,31,41],19:[2,3,13,29,31,41],20:[2,3,5,13,29,31,41]},
 47:{4:[3,5],5:[3,11],6:[3,11,43],7:[3,11,43],8:[3],9:[5,11],10:[11,13,41,43],11:[13,19,41,43],12:[13,19],13:[3,19,37,41,43],14:[3,5,19,41,43],15:[3,17,19,41,43],16:[3],17:[2,3,11,19,37,41,43],18:[2,11,19,31,41,43],19:[2,5,11,31,37,41,43],20:[2,11,29,31,41,43],21:[2,11,29,31,41,43],22:[2,3,29,31,41,43],24:[2]},
 49:{13:[19],19:[2,11,23,41,43,47],23:[2,29,31,41,43,47],24:[2,3]},
 53:{7:[2,5,7,17],8:[2,5,7,17,47],10:[2,5,7,17,23,47],11:[2,5,7,17,23,47],12:[2],16:[5,7],17:[5,7,19,23,41,43,47],19:[5,7,23,41,43,47],20:[5,7],21:[5,7,17,23,41,43,47],22:[5,7],23:[2,5,7,17,41,43,47],24:[2,5,7,17,31,41,47],25:[2,7,17,31,41,43,47]},
 59:{6:[2,7],7:[2,3,11,19],8:[3,11,19],9:[3,11,13,19,53],10:[3],11:[3,5,13,17,19,53],12:[3,5,7,13,17,19,53],13:[2,3,5,7,17,19,53],14:[2,3,5,7,17,19,47,53],15:[2,3,5,7,17,19,23,47,53],16:[3,5,7,17,19,23,47,53],17:[3,5,7,11,19,23,47,53],18:[3,5],19:[3,5,7,11,23,43,47,53],20:[3,5,7,11,23,41,43,47,53],21:[2,3,5,7,11,23,41,43,47,53],22:[2,3,5,7,13,23,41,43,47,53],23:[2,3,5,7,13,19,41,43,47,53],25:[3,7,13,41],26:[3,7,29,41,43,47,53],27:[7,17,19,41,43,47,53],28:[7,41]},
}
for q, d in T2.items():
    for k, ps in d.items():
        for N in ps: DIRECT.append((N, 3, k, q + 1, "LOTW07"))
# semiregular families (orbits of equal length)
def pfac(q):
    for p in range(2, q + 1):
        if q % p == 0: return p
for q in range(7, VMAX + 1):
    if is_pp(q) and q % 4 == 3:
        p0 = pfac(q)
        from math import gcd
        for k in range(3, q - 2):                                  # ASL(1,q), t = 2
            if k % p0 and gcd(k * (k - 1), (q - 1) // 2) == 1:
                DIRECT.append((comb(q - 2, k - 2) // comb(k, 2), 2, k, q, "ASL(1,q)"))
        if q + 1 <= VMAX:
            for k in range(4, q - 2):                              # PSL(2,q), t = 3 (Cusack and Magliveras)
                if gcd(k * (k - 1), p0) == 1 and gcd(k * (k - 1) * (k - 2), (q - 1) // 2) == 1 and gcd(k, (q + 1) // 2) == 1:
                    DIRECT.append((2 * comb(q + 1, k) // (q * (q * q - 1)), 3, k, q + 1, "CuMa99"))
for v in (32, 64):                                                 # AGL(1,2^n), t = 2
    from math import gcd
    for k in range(3, v - 2, 2):
        if gcd(k * (k - 1), v - 1) == 1: DIRECT.append((comb(v - 2, k - 2) // (k * (k - 1)), 2, k, v, "AGL(1,2^n)"))

# orbits of two lengths |G| and |G|/r (LMW01, Theorem 5): affine groups for t = 2, projective groups for t = 3
def divs(n): return [d for d in range(2, n + 1) if n % d == 0]
def is_prime(n): return n > 1 and all(n % d for d in range(2, int(n ** .5) + 1))
for q in range(5, VMAX + 1):
    if not is_pp(q): continue
    p0 = pfac(q)
    for name, m in [("AGL", q - 1)] + ([("ASL", (q - 1) // 2)] if q % 4 == 3 else []):
        order = q * m
        for k in range(3, q - 2):
            if k % p0 == 0: continue
            D = [d for d in divs(m) if k % d == 0 or (k - 1) % d == 0]
            if (not D or (len(D) == 1 and is_prime(D[0]))) and comb(q, k) % order == 0:
                DIRECT.append((comb(q, k) // order, 2, k, q, "LMW01 Thm 5, affine group"))
    if q + 1 <= VMAX:
        for name, m1, m2 in [("PGL", q - 1, q + 1)] + ([("PSL", (q - 1) // 2, (q + 1) // 2)] if q % 4 == 3 else []):
            order = (q + 1) * q * m1
            for k in range(4, q - 2):
                if k % p0 == 0 or (k - 1) % p0 == 0: continue
                D = sorted(set([d for d in divs(m1) if k % d == 0 or (k - 1) % d == 0 or (k - 2) % d == 0] + [d for d in divs(m2) if k % d == 0]))
                if (not D or (len(D) == 1 and is_prime(D[0]) and D[0] > 2)) and comb(q + 1, k) % order == 0:
                    DIRECT.append((comb(q + 1, k) // order, 3, k, q + 1, "LMW01 Thm 5, projective group"))
OURS = [((p - 3) // 4, 3, 4, p, "this paper") for p in (47, 59, 71, 83)] + [(4, 3, 5, 19, "this paper"), (19, 3, 5, 23, "this paper")]
NONEX = [(5,2,3,7,"Cayley"),(5,2,4,7,"complement"),(5,3,4,8,"no LS[5](2,3,7)"),(7,3,4,10,"Kramer and Mesner"),
         (11,4,5,15,"Mendelsohn and Hung"),(13,4,5,17,"Ostergard and Pottonen")]
def axioms(T):
    """General existence theorems (summary table of LOTW07 and the lecture notes KT04)."""
    p = T.p
    for v in range(4, VMAX + 1):
        for k in range(3, v):
            for t in range(2, TMAX + 1):
                if not adm(p, t, k, v) or (p, t, k, v) in NOSET: continue
                kk = min(k, v - k)
                src = None
                if t == 2 and kk == 3: src = "Lu, Teirlinck"
                elif p == 2 and t == 2: src = "Ajoodani-Namini"
                elif p == 2 and t <= 5 and kk <= 15: src = "halving theorems"
                elif p == 2 and t == 6 and kk in (7, 8, 9): src = "halving theorems"
                elif p == 3 and t in (2, 3) and kk <= 80: src = "KT03, LOTW07"
                elif p == 3 and t == 4 and kk <= 10: src = "LOTW07"
                elif p == 5 and t == 2 and kk <= 8: src = "LOTW07"
                elif p == 5 and t == 3 and kk <= 6: src = "LOTW07"
                elif p == 7 and t == 2 and kk <= 6: src = "LMW01"
                elif p == 7 and t == 3 and kk <= 6: src = "LOTW07"
                elif p == 11 and t == 2 and kk <= 10: src = "LMW01"
                elif p in (13, 17, 19, 23) and t == 2 and kk <= 4: src = "LOTW07"
                elif p == 29 and t == 2 and kk <= 28: src = "LOTW07"
                elif t == 3 and kk == 4 and v % 3 == 0 and ((v - 3) // lam_star(3, 4, v)) % p == 0: src = "Teirlinck84"
                if src and k == kk: T.add(t, k, v, src)
NOSET = set()
BAD = set()
def build(with_ours=True):
    tabs = {}
    global NOSET
    NOSET = {(p, t, k, v) for (p, t, k, v, s) in NONEX}
    for p in PRIMES:
        T = Table(p)
        for (q, t, k, v, s) in NONEX:
            if q == p: T.NO[(t, k, v)] = s
        T.close_no()
        NOSET |= {(p,) + x for x in T.NO}
        axioms(T)
        for (N, t, k, v, s) in DIRECT + (OURS if with_ours else []):
            if N % p == 0 and v <= VMAX:
                if not adm(p, t, k, v):
                    BAD.add((N, p, t, k, v, s)); continue
                T.add(t, k, v, s if N == p else f"{s}, LS[{N}]")
        T.close(); T.close_no()
        tabs[p] = T
    return tabs
if __name__ == "__main__":
    import time; t0 = time.time()
    A = build(True); B = build(False)
    for p in PRIMES:
        tot = {t: sum(1 for v in range(t + 2, VMAX + 1) for k in range(t + 1, v // 2 + 1) if adm(p, t, k, v) and v-k>=t) for t in range(2, TMAX + 1)}
        ex = {t: sum(1 for (tt, k, v) in A[p].E if tt == t and 2 * k <= v) for t in range(2, TMAX + 1)}
        new = {t: sum(1 for (tt, k, v) in A[p].E if tt == t and 2 * k <= v and (tt, k, v) not in B[p].E) for t in range(2, TMAX + 1)}
        print(p, " ".join(f"t={t}: {ex[t]}/{tot[t]} (+{new[t]})" for t in range(2, TMAX + 1)))
    print(time.time() - t0)
    for b in sorted(BAD): print('inadmissible base entry dropped:', b)
