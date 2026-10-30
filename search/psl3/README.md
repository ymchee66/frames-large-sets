# Frames with P = PSL(3,q) or kernel containing PSL(3,q)

* `psl3.py` – rank-one frames with H = PSL(3,q) on PG(2,q) ∪ F, q ∈ {2,3,5}, K a point or line stabiliser.
* `n1cert.py` – certifies that each of the 79 open targets in `t_psl3_all.json` fails Proposition 4.1
  (a k-subset that is a union of orbits of the monomial group, which fixes no point and no line, for q ∈ {3,5}).
* `pgammal_cosets.g` – GAP: every nontrivial coset of PSL(3,q) in PΓL(3,q), q ∈ {4,7,8,9}, contains an
  element stabilising k-subsets of PG(2,q) of every size 3 ≤ k ≤ v/2, so frames whose kernel contains
  PSL(3,q) fail Proposition 4.1.
