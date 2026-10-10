import sys, random; sys.path.insert(0,'.')
from v2_words import *
from fractions import Fraction as Fr
def push(w, c, d):
    """apply word (rightmost letter first) to [c,d]; return image if affine on it, else None."""
    for ch in reversed(w):
        g=GEN[ch]; ok=False
        for (x0,y0),(x1,y1) in zip(g,g[1:]):
            if x0<=c and d<=x1:
                sl=(y1-y0)/(x1-x0); c,d=y0+sl*(c-x0), y0+sl*(d-x0); ok=True; break
        if not ok: return None
    return c,d
def check(I,J):
    I,J,Td,Tr=transport(I,J)
    ed,er=exps(Td),exps(Tr)
    w=freered(pos_word(er)+winv(pos_word(ed)))
    return push(w,*addr_iv(I))==addr_iv(P0) and push(w,*addr_iv(J))==addr_iv(Q0)
# cross-check push vs full evaluation on small cases
random.seed(7)
def U(m): return ''.join(random.choice('01') for _ in range(m))
def rnd_pair(m):
    u=U(m); kind=random.choice('SPM')
    if kind=='S':
        v,v2=U(m),U(m)
        while v2==v: v2=U(m)
        return u+'01'+v+'01', u+'01'+v2+'01'
    if kind=='P':
        w,w2=U(m),U(m)
        while w2==w: w2=U(m)
        return w+'01', w2+'01'
    x=U(m)
    while x==u: x=U(m)
    return u+'01'+U(m)+'01', x+'01'
for k in range(30):
    I,J=rnd_pair(3); r=analyze(I,J,check=True); assert r['ok']==check(I,J)==True
print("push-check agrees with full PL evaluation on 30 m=3 cases",flush=True)
# a negative control: wrong-orientation word must fail
I,J,Td,Tr=transport('0001','1101'); w=freered(pos_word(exps(Td))+winv(pos_word(exps(Tr))))
print("negative control (dom*ran^-1) maps correctly?", push(w,*addr_iv(I))==addr_iv(P0), flush=True)
fails=0
for m,ns in ((4,20000),(8,20000),(16,5000),(24,20000),(40,3000),(64,1000)):
    for _ in range(ns):
        I,J=rnd_pair(m); r=analyze(I,J)
        if not (r['F1'] and r['F2a'] and r['F2b'] and check(I,J)): fails+=1; print("FAIL",I,J,r)
    print(f"m={m}: {ns} random pairs: (F1),(F2) and exact transport check; cumulative failures {fails}",flush=True)
