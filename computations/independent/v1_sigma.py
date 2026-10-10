#!/usr/bin/env python3
"""Independent check (mpmath.iv, NOT flint) that sigma(1/32) := inf_{|t-s|>=1/32} ||gamma(t)-gamma(s)|| >= 1/256.
dist^2 = (a_t-a_s)^2 + 2 a_t a_s (1 - sinc(theta_t-theta_s)).
Analytic regions: t<=1 straight; s<=0<1<t uses sinc>=-0.2173 (checked here too); 0<=s,t<=2 linear a>=0;
s>=3 uses monotone sinc; s in[0,3], t>5 uses |sinc|<=1/2.  Remaining compact region s in [0,3], t in [2,5]: BnB."""
from mpmath import iv, mpf
iv.prec = 80
T2 = (iv.mpf(1)/256)**2
def a(t):   # t interval contained in one piece
    lo, hi = t.a, t.b
    if hi <= 2: return t/8
    if lo >= 3: return iv.mpf(5)/16
    assert lo >= 2 and hi <= 3
    return iv.mpf(5)/16 - (3-t)**2/16
def th(t):
    lo, hi = t.a, t.b
    if hi <= 1: return iv.mpf(0)
    if lo >= 2: return t - iv.mpf(3)/2
    assert lo >= 1 and hi <= 2
    return (t-1)**2/2
def pieces(lo, hi, bps):
    pts = [lo]+[b for b in bps if lo < b < hi]+[hi]
    return [iv.mpf([pts[i], pts[i+1]]) for i in range(len(pts)-1)]
def hull(xs):
    return iv.mpf([min(x.a for x in xs), max(x.b for x in xs)])
def A(lo,hi): return hull([a(p) for p in pieces(lo,hi,(2,3))])
def TH(lo,hi): return hull([th(p) for p in pieces(lo,hi,(1,2))])
def oms(q):  # 1 - sinc(q), enclosure
    m = max(abs(q.a), abs(q.b)); mn = 0 if q.a <= 0 <= q.b else min(abs(q.a), abs(q.b))
    if m <= 1:
        return iv.mpf([mpf(mn)**2/6 - mpf(m)**4/120, mpf(m)**2/6])
    return 1 - iv.sin(q)/q
def lb(s0,s1,t0,t1):
    da = A(t0,t1) - A(s0,s1)
    if t1 <= 2: da = iv.mpf([t0-s1, t1-s0])/8
    v = da**2 + 2*A(t0,t1)*A(s0,s1)*oms(TH(t0,t1)-TH(s0,s1))
    return v
def bnb():
    from fractions import Fraction as F
    stack=[(mpf(i)/16, mpf(i+1)/16, mpf(2)+mpf(j)/16, mpf(2)+mpf(j+1)/16) for i in range(48) for j in range(48)]
    n=0; worst=None
    while stack:
        s0,s1,t0,t1 = stack.pop()
        if t1 - s0 < mpf(1)/32: continue          # entire box has t-s<1/32
        # clip: require t-s>=1/32; if box straddles, we still bound over whole box only if lower bound ok
        n+=1
        v = lb(s0,s1,t0,t1)
        if v.a > T2.b:
            worst = v.a if worst is None else min(worst, v.a); continue
        if t0 - s1 < mpf(1)/32 and (s1-s0) < mpf(1)/2**30:
            raise SystemExit(f"stuck {s0} {t0}")
        sm=(s0+s1)/2; tm=(t0+t1)/2
        stack += [(s0,sm,t0,tm),(sm,s1,t0,tm),(s0,sm,tm,t1),(sm,s1,tm,t1)]
    return n, worst
if __name__=="__main__":
    # sinc lower bound on [1, 1/0.2173]
    st=[(mpf(1), mpf(1)/mpf('0.2173'))]; k=0
    while st:
        l,h=st.pop(); k+=1; q=iv.mpf([l,h]); v=iv.sin(q)/q
        if v.a > -0.2173: continue
        m=(l+h)/2; st+=[(l,m),(m,h)]
    print("sinc>=-0.2173 on [1,4.602]: boxes",k)
    print("R4: 2(5/16)^2(1-sinc(1/32)) lower", (2*(iv.mpf(5)/16)**2*oms(iv.mpf(1)/32)).a, " target", T2.b)
    n,w = bnb(); print("R6 boxes", n, "min lower bound dist^2", w, " target", T2.b, " ratio", w/T2.b)
