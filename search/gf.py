"""Small finite fields GF(p^n) as integer-coded polynomials; enough for q<=343."""
from itertools import product
def is_prime(n): return n>1 and all(n%d for d in range(2,int(n**.5)+1))
def factor_pp(q):
    for p in range(2,q+1):
        if is_prime(p) and q%p==0:
            n=0; m=q
            while m%p==0: m//=p; n+=1
            return (p,n) if m==1 else None
class GF:
    def __init__(s,q):
        s.q=q; s.p,s.n=factor_pp(q); p,n=s.p,s.n
        if n==1:
            s.add=[[(a+b)%p for b in range(p)] for a in range(p)]
            s.mulT=[[(a*b)%p for b in range(p)] for a in range(p)]
        else:
            def dig(a): return [(a//p**i)%p for i in range(n)]
            def und(d): return sum(x*p**i for i,x in enumerate(d))
            # find irreducible monic poly of degree n by brute force (no roots is not enough for n>3; test by field property)
            for coeffs in product(range(p),repeat=n):
                f=list(coeffs)+[1]
                def red(c):
                    c=c[:]
                    for i in range(len(c)-1,n-1,-1):
                        t=c[i]
                        if t:
                            for j in range(n+1): c[i-n+j]=(c[i-n+j]-t*f[j])%p
                    return c[:n]
                def pm(a,b):
                    A=dig(a);B=dig(b);c=[0]*(2*n-1)
                    for i in range(n):
                        if A[i]:
                            for j in range(n): c[i+j]=(c[i+j]+A[i]*B[j])%p
                    return und(red(c))
                # primitive element test: x (=p) has order q-1
                x=p; y=x; o=1
                while y!=1 and o<q: y=pm(y,x); o+=1
                if o==q-1:
                    s.mulT=[[pm(a,b) for b in range(q)] for a in range(q)]
                    s.add=[[und([(u+w)%p for u,w in zip(dig(a),dig(b))]) for b in range(q)] for a in range(q)]
                    break
            else: raise ValueError
        s.neg=[next(b for b in range(q) if s.add[a][b]==0) for a in range(q)]
        s.inv=[None]+[next(b for b in range(1,q) if s.mulT[a][b]==1) for a in range(1,q)]
        s.squares={s.mulT[a][a] for a in range(1,q)}
    def frob(s,a):
        r=1
        for _ in range(s.p): r=s.mulT[r][a]
        return r
