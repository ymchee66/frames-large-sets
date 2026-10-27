def gf(p,m):
    import itertools as it
    for c in it.product(range(p),repeat=m):
        f=list(c)+[1]   # monic x^m + ...
        q=p**m
        def mul(a,b):
            A=[(a//p**i)%p for i in range(m)]; B=[(b//p**i)%p for i in range(m)]
            P=[0]*(2*m-1)
            for i in range(m):
                for j in range(m): P[i+j]=(P[i+j]+A[i]*B[j])%p
            for d in range(2*m-2,m-1,-1):
                co=P[d]
                if co:
                    for i in range(m+1): P[d-m+i]=(P[d-m+i]-co*f[i])%p
            return sum(P[i]*p**i for i in range(m))
        def add(a,b): return sum((((a//p**i)+(b//p**i))%p)*p**i for i in range(m))
        # primitive element test: x (=p) generates
        x=p if m>1 else 2; y=1; seen=set()
        for _ in range(q-1): y=mul(y,x); seen.add(y)
        if len(seen)==q-1 and 0 not in seen: return mul,add,x
    raise ValueError
