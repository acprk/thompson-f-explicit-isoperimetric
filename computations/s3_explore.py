#!/usr/bin/env python3
"""NON-RIGOROUS exploration (Monte Carlo sampling). Average transport-word lengths for two
alternative couplings ("push"/"pull") of descendant/parent pairs. Nothing in the note depends on it."""
import sys, random; sys.path.insert(0,'.')
from thompson import *
def P(u): return u+"01"
def Q(u,v): return u+"01"+v+"01"
def couple(dom, ran):
    (c,d),(a,b)=sorted(dom,key=lo_end),sorted(ran,key=lo_end)
    # inner product symmetric: lower->lower
    res=suffix_transport(c,d,a,b)
    return res, (c,d,a,b)
def stats(m, samples, seed=0, verify=0):
    random.seed(seed); U=[format(i,'0%db'%m) for i in range(2**m)]
    acc={"push":[], "pull":[]}
    nver=0
    for _ in range(samples):
        u=random.choice(U); w=random.choice(U); x=random.choice(U)
        if x in (u,w): continue
        for fam in ("push","pull"):
            if fam=="push": dom=[Q(u,w),P(x)]; ran=[P(w),P(x)]
            else: dom=[Q(u,w),P(x)]; ran=[Q(u,w),Q(u,x)]
            res,(c,d,a,b)=couple(dom,ran)
            (Td,Tr),N=res
            wd=nf_word(Td,Tr); n0,n1=letter_counts(wd)
            if nver<verify:
                assert check_transport(wd,c,d,a,b); nver+=1
            acc[fam].append((N,len(wd),n0,n1))
    for fam,v in acc.items():
        k=len(v)
        print(f"m={m} r={m+2} {fam}: avg N {sum(t[0] for t in v)/k:.2f} max N {max(t[0] for t in v)}  avg|w| {sum(t[1] for t in v)/k:.2f} avg n0 {sum(t[2] for t in v)/k:.2f} avg n1 {sum(t[3] for t in v)/k:.2f} max|w| {max(t[1] for t in v)}")
for m in (4,8,12,20):
    stats(m, 300 if m<20 else 150, verify=40 if m<=8 else 5)
