import sys, itertools; sys.path.insert(0,'.')
from fractions import Fraction as Fr
from v2_words import analyze
def closed(m):
    D=2**m; r=m+2; lam=Fr(D-m-1,D-1)
    Nbar=Fr(11*r,4)-1-lam
    def lzge(k): return D//2**k if k<=m else (1 if k==m+1 else 0)
    ES=sum(Fr(lzge(k),D) for k in range(1,m+2))
    def maxge(k):
        nb=D-lzge(k); return Fr(D*(D-1)-nb*(nb-1),D*(D-1))
    EP=sum(maxge(k) for k in range(1,m+2))
    e0=(ES+EP+2*EP)/4
    return Nbar,e0,max(2*Nbar+e0-4,2*Nbar-2)
def exhaustive(m):
    U=[''.join(t) for t in itertools.product('01',repeat=m)]; D=len(U)
    cache={}
    def A(I,J):
        k=(I,J)
        if k not in cache: cache[k]=analyze(I,J)
        return cache[k]
    s={'N':0,'n0':0,'n1':0,'e0':0,'W':0}
    Pstats={'N':0,'n0':0,'n1':0,'e0':0,'W':0}
    for w,w2 in itertools.permutations(U,2):
        r=analyze(w+'01',w2+'01')
        for key in ('N','n0','n1'): Pstats[key]+=r[key]
        Pstats['e0']+=r['e0d']; Pstats['W']+=1
        assert r['F1'] and r['F2a'] and r['F2b']
    for u in U:
        for v,v2 in itertools.permutations(U,2):
            r=analyze(u+'01'+v+'01',u+'01'+v2+'01'); wt=1
            assert r['F1'] and r['F2a'] and r['F2b']
            for key in ('N','n0','n1'): s[key]+=wt*r[key]
            s['e0']+=wt*r['e0d']; s['W']+=wt
        for key in Pstats: s[key]+=Pstats[key]
        for v in U:
            for x in U:
                if x==u: continue
                r=analyze(u+'01'+v+'01',x+'01'); wt=2
                assert r['F1'] and r['F2a'] and r['F2b']
                for key in ('N','n0','n1'): s[key]+=wt*r[key]
                s['e0']+=wt*r['e0d']; s['W']+=wt
    W=s['W']
    return {k:Fr(v,W) for k,v in s.items() if k!='W'}, W
if __name__=="__main__":
    for m in range(1,(int(sys.argv[1]) if len(sys.argv)>1 else 4)+1):
        av,W=exhaustive(m); Nb,e0,M=closed(m)
        print(f"m={m} W={W}=4D^2(D-1)? {W==4*4**m*(2**m-1)} avgN={float(av['N']):.4f} closed={float(Nb):.4f} eq={av['N']==Nb} "
              f"avg e0={float(av['e0']):.4f} closed={float(e0):.4f} eq={av['e0']==e0} avg n0={float(av['n0']):.3f} avg n1={float(av['n1']):.3f} Mbar={float(M):.3f} "
              f"ok={max(av['n0'],av['n1'])<=M}", flush=True)
