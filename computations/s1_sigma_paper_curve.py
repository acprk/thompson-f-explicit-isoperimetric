#!/usr/bin/env python3
"""
Make the ONE non-explicit constant of the paper's displacement map explicit.

In the Lean formalization (AnalyticConstruction.lean L22-30) the tube radius is
    rho = min(1/1024, sigma/8),  sigma = Classical.choose (gamma_uniform_separation (1/32))
(TubeRegularCurve.lean L46), i.e. sigma is ANY positive lower bound for
    inf_{|t-s| >= 1/32} || gamma(t) - gamma(s) ||_{L^2([0,1];R^2)}.
The straight part (t,s <= 1, gamma = (t/8) v0) gives sigma <= 1/256.  We prove sigma >= 1/256
(hence sigma = 1/256 exactly) rigorously:

  ||gamma(t)-gamma(s)||^2 = a(t)^2 + a(s)^2 - 2 a(t) a(s) sinc(theta(t)-theta(s))   (auxiliary.tex L94)
                          = (a(t)-a(s))^2 + 2 a(t) a(s) (1 - sinc(.)).

Regions (s<t, t-s>=1/32):
 R1 t<=1                : theta=0, distance=(t-s)/8 >= 1/256 (equality case).
 R2 s<=0<1<t            : a(s)<=0<1/8<=a(t); sinc >= -0.2173  => dist^2 >= (1-0.2173)(a_t^2+a_s^2) >= 0.78/64.
 R3 0<=s, t<=2          : a linear (slope 1/8) and a_t a_s >= 0  => dist >= (t-s)/8.
 R4 s>=3                : dist^2 = 2(5/16)^2 (1 - sinc(t-s)) >= 2(5/16)^2(1-sinc(1/32))  (sinc decreasing on [0,pi],
                          |sinc|<=1/pi < sinc(1/32) beyond).
 R5 0<=s<=3, t>5        : theta_t-theta_s >= t-3 >= 2, |sinc|<=1/2 => dist^2 >= a_t^2 (1-1/4) >= (1/4)^2*3/4.
 R6 0<=s<=3, 2<t<=5     : compact; rigorous branch-and-bound with Arb balls in variables (s, h=t-s).
All arithmetic below in python-flint arb (outward rounding).
"""
from flint import arb as _arb, ctx

def arb(m, r=None):
    # inflate radii so float box endpoints are always enclosed (outward safety margin)
    if r is None:
        return _arb(m)
    return _arb(m, abs(r) * (1 + 1e-9) + 1e-300)
ctx.prec = 80

ONE = arb(1)
def _split(t_lo, t_hi, bps):
    pts = [t_lo] + [b for b in bps if t_lo < b < t_hi] + [t_hi]
    return list(zip(pts[:-1], pts[1:]))

def a_iv(t_lo, t_hi):
    parts = [_a_iv(l, h) for l, h in _split(t_lo, t_hi, (2, 3))]
    out = parts[0]
    for p in parts[1:]:
        out = out.union(p)
    return out

def theta_iv(t_lo, t_hi):
    parts = [_theta_iv(l, h) for l, h in _split(t_lo, t_hi, (1, 2))]
    out = parts[0]
    for p in parts[1:]:
        out = out.union(p)
    return out

def _a_iv(t_lo, t_hi):
    """enclosure of a on [t_lo,t_hi] (box must not straddle 2 or 3)."""
    t = arb((t_lo + t_hi) / 2, (t_hi - t_lo) / 2)
    if t_hi <= 2:
        return t / 8
    if t_lo >= 3:
        return arb(5) / 16
    assert t_lo >= 2 and t_hi <= 3
    return arb(5) / 16 - (3 - t) ** 2 / 16

def _theta_iv(t_lo, t_hi):
    t = arb((t_lo + t_hi) / 2, (t_hi - t_lo) / 2)
    if t_hi <= 1:
        return arb(0)
    if t_lo >= 2:
        return t - arb(3) / 2
    assert t_lo >= 1 and t_hi <= 2
    return (t - 1) ** 2 / 2

def one_minus_sinc(q):
    """enclosure of 1 - sin(q)/q for an arb q (handles q near 0 via series bound)."""
    # 1 - sinc(q) = q^2/6 - q^4/120 + ...; for |q| <= 1 use q^2/6 - q^4/120 <= . <= q^2/6
    mid = abs(q.mid()); rad = q.rad()
    qa = abs(q)
    if float(mid) + float(rad) <= 1.0:
        lo = qa.lower() ** 2 / 6 - qa.upper() ** 4 / 120
        hi = qa.upper() ** 2 / 6
        return arb((lo + hi) / 2, (hi - lo) / 2)
    return ONE - q.sin() / q

TARGET = (arb(1) / 256) ** 2

def lower_dist2(s_lo, s_hi, h_lo, h_hi):
    t_lo, t_hi = s_lo + h_lo, s_hi + h_hi
    # split t-range at breakpoints if needed -> caller ensures alignment by bisection
    as_ = a_iv(s_lo, s_hi); at = a_iv(t_lo, t_hi)
    # a_t - a_s with dependency reduction: if both in linear region use h/8
    if t_hi <= 2 and s_hi <= 2:
        da = arb((h_lo + h_hi) / 2, (h_hi - h_lo) / 2) / 8
    else:
        da = at - as_
    dth = theta_iv(t_lo, t_hi) - theta_iv(s_lo, s_hi)
    val = da ** 2 + 2 * at * as_ * one_minus_sinc(dth)
    return val

def bnb():
    from fractions import Fraction as Fr
    # region R6 in (s,h): s in [0,3], h in [1/32, 5-s], t=s+h in (2,5]
    stack = []
    N = 48
    for i in range(N):
        for j in range(4 * N):
            s_lo, s_hi = 3 * i / N, 3 * (i + 1) / N
            h_lo, h_hi = 1 / 32 + 5 * j / (4 * N), 1 / 32 + 5 * (j + 1) / (4 * N)
            stack.append((s_lo, s_hi, h_lo, h_hi, 0))
    boxes = 0; worst = None
    while stack:
        s_lo, s_hi, h_lo, h_hi, dep = stack.pop()
        t_lo, t_hi = s_lo + h_lo, s_hi + h_hi
        if t_hi <= 2 or t_lo > 5:
            continue  # outside R6 (R3 / R5)
        boxes += 1
        v = lower_dist2(s_lo, s_hi, h_lo, h_hi)
        if v is not None and v > TARGET:
            lo = v.lower()
            if worst is None or lo < worst:
                worst = lo
            continue
        if dep > 40:
            raise RuntimeError(f"bnb failed near s=[{s_lo},{s_hi}] h=[{h_lo},{h_hi}] v={v}")
        sm, hm = (s_lo + s_hi) / 2, (h_lo + h_hi) / 2
        # bisect the wider (in t-terms) coordinate
        stack += [(s_lo, sm, h_lo, hm, dep + 1), (sm, s_hi, h_lo, hm, dep + 1),
                  (s_lo, sm, hm, h_hi, dep + 1), (sm, s_hi, hm, h_hi, dep + 1)]
    return boxes, worst

def check_sinc_min():
    # sinc(q) >= -0.2173 on [0, 1/0.2173] by bnb; beyond, |sinc| <= 1/q <= 0.2173.
    lim = 1 / 0.2173
    stack = [(1.0, lim)]
    n = 0
    while stack:
        lo, hi = stack.pop(); n += 1
        q = arb((lo + hi) / 2, (hi - lo) / 2)
        v = q.sin() / q
        if v > arb(-0.2173):
            continue
        m = (lo + hi) / 2; stack += [(lo, m), (m, hi)]
    return n  # [0,1]: sinc>=sin(1)>0

if __name__ == "__main__":
    n = check_sinc_min()
    print("R2 lemma: sinc >= -0.2173 verified, boxes", n)
    r2 = (1 - arb(0.2173)) / 64
    print("R2 bound dist^2 >=", r2, "> target", bool(r2 > TARGET))
    r4 = 2 * (arb(5) / 16) ** 2 * one_minus_sinc(arb(1) / 32)
    print("R4 bound dist^2 >=", r4.lower(), "> target", bool(r4 > TARGET))
    r5 = (arb(1) / 4) ** 2 * arb(3) / 4
    print("R5 bound dist^2 >=", r5, bool(r5 > TARGET))
    boxes, worst = bnb()
    print("R6 branch-and-bound boxes:", boxes, " min lower bound of dist^2:", worst,
          " target (1/256)^2 =", TARGET.mid())
    print("=> sigma(1/32) = 1/256 exactly (upper bound from straight segment t,s<=1).")
    # paper/Lean constants made explicit
    from fractions import Fraction as Fr
    sigma = Fr(1, 256)
    rho = min(Fr(1, 1024), sigma / 8)
    LipG = max(Fr(1), Fr(30) + Fr(7) / (4 * rho) + 3808 * rho, Fr(5) + Fr(7) / (8 * rho))
    L = 2 * LipG / (rho / 4)
    print("rho =", rho, " Lip(G) =", LipG, "=", float(LipG), " L_paper =", L, "=", float(L))
