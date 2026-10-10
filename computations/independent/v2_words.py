#!/usr/bin/env python3
"""Independent re-implementation: Thompson F as exact PL maps (integer numerators over 2^K),
CFP normal form from tree pairs, {x0,x1} words, transport check, and (F1)/(F2) counts.
Convention hg = h o g ; word s1...sk evaluates to s1 o ... o sk."""
import sys, random, itertools
from fractions import Fraction as Fr

# ---- PL maps as sorted breakpoint lists of Fractions ----
def ev(bp, x):
    for (x0,y0),(x1,y1) in zip(bp, bp[1:]):
        if x0 <= x <= x1:
            return y0 + (y1-y0)*(x-x0)/(x1-x0)
    raise ValueError
def inv(bp): return [(y,x) for x,y in bp]
def comp(h, g):  # h o g
    xs = sorted(set(x for x,_ in g) | set(ev(inv(g), x) for x,_ in h))
    return simplify([(x, ev(h, ev(g, x))) for x in xs])
def simplify(bp):
    out=[bp[0]]
    for i in range(1,len(bp)-1):
        (x0,y0),(x1,y1),(x2,y2) = out[-1], bp[i], bp[i+1]
        if (y1-y0)*(x2-x1) != (y2-y1)*(x1-x0): out.append(bp[i])
    out.append(bp[-1]); return out
H=Fr(1,2)
X0=[(Fr(0),Fr(0)),(Fr(1,2),Fr(1,4)),(Fr(3,4),Fr(1,2)),(Fr(1),Fr(1))]
X1=[(Fr(0),Fr(0)),(Fr(1,2),Fr(1,2)),(Fr(3,4),Fr(5,8)),(Fr(7,8),Fr(3,4)),(Fr(1),Fr(1))]
GEN={'a':X0,'A':inv(X0),'b':X1,'B':inv(X1)}
def evword(w):
    g=[(Fr(0),Fr(0)),(Fr(1),Fr(1))]
    for c in reversed(w): g = comp(GEN[c], g)
    return g
def addr_iv(w):
    k=int(w,2) if w else 0; n=len(w); return Fr(k,2**n), Fr(k+1,2**n)
def affine_onto(g, I, J):
    (a,b),(c,d)=addr_iv(I),addr_iv(J)
    if ev(g,a)!=c or ev(g,b)!=d: return False
    return all(not (a<x<b) for x,_ in g)
# ---- trees: set of leaf addresses ----
def mintree(addrs):
    internal=set()
    for w in addrs:
        for i in range(len(w)): internal.add(w[:i])
    return sorted({v+c for v in internal for c in '01' if v+c not in internal}, key=lambda w: addr_iv(w)[0])
def exps(leaves):
    """CFP leaf exponents (Cannon-Floyd-Parry §2): length of maximal left-edge path from the leaf
    not reaching the right side of the tree."""
    out=[]
    for w in leaves:
        e=len(w)-len(w.rstrip('0')); top=w[:len(w)-e]
        if e>0 and set(top)<= {'1'}: e-=1   # top vertex on right spine: that edge touches right side
        out.append(e)
    return out
def pos_word(e):
    """x_0^{e0} x_1^{e1} ... x_n^{en} with x_k = x0^{-(k-1)} x1 x0^{k-1} (k>=1), telescoped."""
    n=max([k for k in range(1,len(e)) if e[k]>0], default=0)
    w='a'*e[0]
    for k in range(1,n+1):
        w+='b'*e[k]
        if k<n: w+='A'
    if n>=1: w+='a'*(n-1)
    return w
def winv(w): return ''.join(c.swapcase() for c in reversed(w))
def freered(w):
    s=[]
    for c in w:
        if s and s[-1]==c.swapcase(): s.pop()
        else: s.append(c)
    return ''.join(s)
def counts(w): return sum(c in 'aA' for c in w), sum(c in 'bB' for c in w)
REF=['00','01','100','101','11']; P0,Q0='01','101'
def pad_last(gap,k):
    gap=list(gap)
    while len(gap)<k:
        w=gap.pop(); gap+= [w+'0',w+'1']
    return gap
def transport(I,J):
    if addr_iv(I)[0]>addr_iv(J)[0]: I,J=J,I
    Td=mintree([I,J]); i,j=Td.index(I),Td.index(J)
    gd=[Td[:i],Td[i+1:j],Td[j+1:]]
    assert all(len(g)>=1 for g in gd)
    gr=[['00'],['100'],['11']]
    pr=[pad_last(y,len(x)) for x,y in zip(gd,gr)]
    Tr=pr[0]+[P0]+pr[1]+[Q0]+pr[2]
    assert len(Tr)==len(Td)
    return I,J,Td,Tr
def lcp(x,y):
    n=0
    while n<min(len(x),len(y)) and x[n]==y[n]: n+=1
    return n
def analyze(I,J,check=False, which="ran*dom^-1"):
    I,J,Td,Tr=transport(I,J)
    N=len(Td)-1
    F1 = (N == len(I)+len(J)-lcp(I,J)-1)
    ed,er=exps(Td),exps(Tr)
    w_un = pos_word(er)+winv(pos_word(ed))
    w=freered(w_un)
    n0,n1=counts(w)
    e0d,e0r=ed[0],er[0]; sd,sr=len(Td[-1]),len(Tr[-1])
    F2a = n1 <= (N-sr-e0r)+(N-sd-e0d)
    F2b = n0 <= e0r+e0d+2*N-6
    ok=None
    if check:
        g=evword(w); ok = affine_onto(g,I,P0) and affine_onto(g,J,Q0)
    return dict(N=N,n0=n0,n1=n1,e0d=e0d,e0r=e0r,F1=F1,F2a=F2a,F2b=F2b,ok=ok,lz=len(I)-len(I.lstrip('0')))
