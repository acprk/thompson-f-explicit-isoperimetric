import sys, random, itertools; sys.path.insert(0,'.')
from v2_words import *
# sanity of generators / normal form on known trees
assert evword('a')==X0
# x2 = x0^-1 x1 x0
random.seed(1)
def rnd_pair(m):
    u=''.join(random.choice('01') for _ in range(m))
    kind=random.choice('SPM')
    U=lambda: ''.join(random.choice('01') for _ in range(m))
    if kind=='S':
        v,v2=U(),U()
        while v2==v: v2=U()
        return u+'01'+v+'01', u+'01'+v2+'01'
    if kind=='P':
        w,w2=U(),U()
        while w2==w: w2=U()
        return w+'01', w2+'01'
    x=U()
    while x==u: x=U()
    return u+'01'+U()+'01', x+'01'
bad=0
for m,ns,nchk in ((5,3000,3000),(8,2000,500),(12,1000,100),(24,1000,40),(40,300,10),(60,200,5)):
    tot=0
    for k in range(ns):
        I,J=rnd_pair(m)
        r=analyze(I,J,check=(k<nchk))
        if not (r['F1'] and r['F2a'] and r['F2b'] and r['ok'] in (None,True)):
            bad+=1; print("FAIL",m,I,J,r)
    print(f"m={m}: {ns} random pairs, {min(ns,nchk)} transports evaluated exactly; failures so far {bad}")
