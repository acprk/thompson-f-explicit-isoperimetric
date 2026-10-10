#!/usr/bin/env python3
"""
Certificate for "An explicit isoperimetric constant for Thompson's group F".

Everything that enters a stated constant is computed either in exact rational arithmetic
(fractions.Fraction) or in Arb ball arithmetic (python-flint, outward rounding).  Floats are
used only for printing.

Parts
  [A] sigma(1/32) = 1/256 for the curve gamma of the OpenAI preprint (App. A), rigorous:
      analytic regions + Arb branch-and-bound on the compact region.
  [B] constant chain rho, Lip(G), L (old 2/eps and Hilbert 1/eps normalisation), kappa.
  [C] word-length bookkeeping for the {x0,x1} conversion: closed forms Nbar(m), e0bar(m), Mbar(m),
      plus an exhaustive unit test (m <= 4) of (F1), (F2), the averages and the transports.
  [D] the main bounds h_E >= G(kappa) = max_m (kappa - 4/2^m) / (4 (1-2^-m) Mbar(m)),
      unconditional (1/eps and 2/eps) and conditional g(k) for a k-Lipschitz retraction B -> S.
  [E] framework ceilings, loss ledger, spectral corollary.

Run:  nice -n 10 python3 cert.py          (about 15 s, single core)
"""
import itertools, math, sys
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx

ctx.prec = 128

def A(x):
    """exact rational -> arb (exact when x is dyadic, otherwise a tight ball)."""
    x = Fr(x)
    return arb(fmpq(x.numerator, x.denominator))

def dec(x, digits=6):
    """decimal string of an exact rational via arb."""
    return A(x).str(digits, radius=False)

# ---------------------------------------------------------------------------------------------
# [A] sigma(1/32) = 1/256
# ---------------------------------------------------------------------------------------------
# a(t) = t/8 (t<=2), 5/16-(3-t)^2/16 (2<t<3), 5/16 (t>=3);  theta(t) = 0 (t<=1), (t-1)^2/2 (1<t<2),
# t-3/2 (t>=2);  ||gamma(t)-gamma(s)||^2 = a_t^2 + a_s^2 - 2 a_t a_s sinc(theta_t - theta_s).

def hull(lo, hi):
    return A(lo).union(A(hi))

def pieces(lo, hi, bps):
    pts = [lo] + [b for b in bps if lo < b < hi] + [hi]
    return list(zip(pts[:-1], pts[1:]))

def a_enc(lo, hi):
    out = None
    for l, h in pieces(lo, hi, (Fr(2), Fr(3))):
        T = hull(l, h)
        if h <= 2:
            v = T / 8
        elif l >= 3:
            v = A(Fr(5, 16))
        else:
            v = A(Fr(5, 16)) - (3 - T) ** 2 / 16
        out = v if out is None else out.union(v)
    return out

def th_enc(lo, hi):
    out = None
    for l, h in pieces(lo, hi, (Fr(1), Fr(2))):
        T = hull(l, h)
        if h <= 1:
            v = A(0)
        elif l >= 2:
            v = T - A(Fr(3, 2))
        else:
            v = (T - 1) ** 2 / 2
        out = v if out is None else out.union(v)
    return out

def one_minus_sinc(q):
    """enclosure of 1 - sin(q)/q.  For |q| <= 1 use q^2/6 - q^4/120 <= 1 - sinc q <= q^2/6
    (alternating series with decreasing terms); the bounds use explicit |q|-range endpoints."""
    lo, hi = q.lower(), q.upper()          # arb endpoints (exact balls of radius 0)
    if hi <= 1 and lo >= -1:
        amax = max(abs(lo), abs(hi))
        amin = A(0) if (lo <= 0 and hi >= 0) else min(abs(lo), abs(hi))
        L = amin ** 2 / 6 - amax ** 4 / 120
        U = amax ** 2 / 6
        return L.union(U)
    return 1 - q.sin() / q

TARGET = A(Fr(1, 256)) ** 2

def lower_dist2(s0, s1, h0, h1):
    t0, t1 = s0 + h0, s1 + h1
    as_, at = a_enc(s0, s1), a_enc(t0, t1)
    if t1 <= 2:                          # both in the linear part: a_t - a_s = h/8 exactly
        da = hull(h0, h1) / 8
    else:
        da = at - as_
    q = th_enc(t0, t1) - th_enc(s0, s1)
    return da ** 2 + 2 * at * as_ * one_minus_sinc(q)

def part_A():
    print("[A] sigma(1/32) = 1/256")
    # A1: sinc(q) >= -0.2173 for all q (q in [0,1]: sinc >= sin 1 > 0; q >= 1/0.2173: |sinc| <= 1/q)
    c = Fr(2173, 10000)
    stack, n = [(Fr(1), 1 / c)], 0
    while stack:
        lo, hi = stack.pop(); n += 1
        q = hull(lo, hi)
        if q.sin() / q > -A(c):
            continue
        assert hi - lo > Fr(1, 2 ** 40)
        m = (lo + hi) / 2; stack += [(lo, m), (m, hi)]
    print(f"  A1 sinc >= -0.2173 on [1, 1/0.2173]: {n} boxes")
    # A2 (s<=0<1<t): dist^2 >= (1-0.2173) a_t^2 >= (1-0.2173)/64
    r2 = (1 - A(c)) / 64
    # A4 (s>=3): dist^2 = 2(5/16)^2(1-sinc(t-s)) >= 2(5/16)^2(1-sinc(1/32))
    r4 = 2 * A(Fr(5, 16)) ** 2 * one_minus_sinc(A(Fr(1, 32)))
    # A5 (0<=s<=3, t>5): theta_t-theta_s >= 2, |sinc|<=1/2, a_t = 5/16 >= 1/4: dist^2 >= (3/4) a_t^2
    r5 = A(Fr(3, 4)) * A(Fr(1, 4)) ** 2
    for name, r in (("A2", r2), ("A4", r4), ("A5", r5)):
        assert r > TARGET
        print(f"  {name}: dist^2 >= {r.lower().str(6)}  > (1/256)^2 = {TARGET.str(6)}  ratio {(r / TARGET).lower().str(4)}")
    # A6: compact region 0<=s<=3, 2<t<=5, coordinates (s, h=t-s), h >= 1/32
    N = 48
    stack = []
    for i in range(N):
        for j in range(4 * N):
            s0, s1 = Fr(3 * i, N), Fr(3 * (i + 1), N)
            h0 = Fr(1, 32) + Fr(5 * j, 4 * N); h1 = Fr(1, 32) + Fr(5 * (j + 1), 4 * N)
            stack.append((s0, s1, h0, h1, 0))
    boxes, worst = 0, None
    while stack:
        s0, s1, h0, h1, d = stack.pop()
        if s1 + h1 <= 2 or s0 + h0 > 5:
            continue                      # entirely in the analytic regions A3 / A5
        boxes += 1
        v = lower_dist2(s0, s1, h0, h1)
        if v > TARGET:
            lo = v.lower()
            worst = lo if worst is None or lo < worst else worst
            continue
        if d > 40:
            raise RuntimeError(f"BnB stuck at s=[{s0},{s1}] h=[{h0},{h1}]")
        sm, hm = (s0 + s1) / 2, (h0 + h1) / 2
        stack += [(s0, sm, h0, hm, d + 1), (sm, s1, h0, hm, d + 1),
                  (s0, sm, hm, h1, d + 1), (sm, s1, hm, h1, d + 1)]
    print(f"  A6 branch-and-bound: {boxes} boxes certified, min certified lower bound of dist^2 = "
          f"{worst.str(6)}  (ratio to target {(worst / TARGET).str(4)})")
    print("  A1-A6 + straight segment (t,s<=1: dist=(t-s)/8, equality at t-s=1/32)  =>  sigma(1/32) = 1/256")
    return boxes

# ---------------------------------------------------------------------------------------------
# [B] constant chain
# ---------------------------------------------------------------------------------------------
SIGMA = Fr(1, 256)
RHO = min(Fr(1, 1024), SIGMA / 8)                                   # TubeRegularCurve.lean: rho
LIPG = max(Fr(1), 30 + Fr(7) / (4 * RHO) + 3808 * RHO, 5 + Fr(7) / (8 * RHO))   # AnalyticCorrection.lean
EPS = RHO / 4                                                       # ||G|| >= rho/4
L_OLD = 2 * LIPG / EPS                                              # AnalyticConstruction.lean (2/eps)
L_NEW = LIPG / EPS                                                  # Hilbert: normalisation is 1/eps-Lipschitz
DELTA = Fr(1, 2)

def part_B():
    print("[B] constant chain")
    print(f"  sigma = {SIGMA}, rho = min(1/1024, sigma/8) = {RHO}")
    print(f"  Lip(G) = max(1, 30+7/(4rho)+3808rho, 5+7/(8rho)) = {LIPG} = {dec(LIPG, 10)}")
    print(f"  eps = rho/4 = {EPS}")
    print(f"  L (2/eps, as in Lean)      = {L_OLD}")
    print(f"  L (1/eps, Hilbert)         = {L_NEW}")
    print(f"  kappa = (delta/L)^2: old {dec((DELTA / L_OLD) ** 2)}, new {dec((DELTA / L_NEW) ** 2)}")
    # numerical sanity of ||x/|x| - y/|y||| <= ||x-y||/sqrt(|x||y|) (identity: difference = (|x|-|y|)^2)
    import random
    random.seed(1)
    worst = 0.0
    for _ in range(20000):
        dmn = random.choice((2, 3, 10))
        x = [random.gauss(0, 1) * random.choice((1e-3, 1, 5)) for _ in range(dmn)]
        y = [random.gauss(0, 1) * random.choice((1e-3, 1, 5)) for _ in range(dmn)]
        nx, ny = math.sqrt(sum(t * t for t in x)), math.sqrt(sum(t * t for t in y))
        lhs = math.sqrt(sum((a / nx - b / ny) ** 2 for a, b in zip(x, y)))
        rhs = math.sqrt(sum((a - b) ** 2 for a, b in zip(x, y))) / math.sqrt(nx * ny)
        worst = max(worst, lhs / rhs if rhs > 0 else 0)
    print(f"  random test of the Hilbert normalisation inequality: max lhs/rhs = {worst:.12f} (must be <= 1)")
    assert worst <= 1 + 1e-12

# ---------------------------------------------------------------------------------------------
# [C] word lengths
# ---------------------------------------------------------------------------------------------
def closed_form_sharp(m):
    """(F2') version: Mbar'(m) = max(2Nbar - 2gbar + e0bar - 2, 2Nbar - gbar - 3 - e0bar),
    gbar = weighted average of #zeros of the upper interval = (3m+4)/4."""
    Nb, e0b, _ = closed_form(m)
    gb = Fr(3 * m + 4, 4)
    return Nb, e0b, gb, max(2 * Nb - 2 * gb + e0b - 2, 2 * Nb - gb - 3 - e0b)

def closed_form(m):
    """exact averages over the weighted family (sibling 1, parent 1, mixed 2) of N and e0_dom,
    and Mbar(m) = max(2 Nbar + e0bar - 4, 2 Nbar - 2)."""
    D = 2 ** m; r = m + 2
    lam = Fr(D - m - 1, D - 1)                          # E[lcp(v,v') | v != v'], v,v' in {0,1}^m
    Nbar = Fr(11 * r, 4) - 1 - lam
    def lz_ge(k):                                       # #{u : lz(u01) >= k}
        return D // 2 ** k if k <= m else (1 if k == m + 1 else 0)
    ES = sum(Fr(lz_ge(k), D) for k in range(1, m + 2))  # = 1
    def max_ge(k):
        nb = D - lz_ge(k)
        return Fr(D * (D - 1) - nb * (nb - 1), D * (D - 1))
    EP = sum(max_ge(k) for k in range(1, m + 2))
    e0bar = (ES + 3 * EP) / 4
    return Nbar, e0bar, max(2 * Nbar + e0bar - 4, 2 * Nbar - 2)

def iv(w):
    k = int(w, 2) if w else 0
    return Fr(k, 2 ** len(w)), Fr(k + 1, 2 ** len(w))

def mintree(addrs):
    internal = {w[:i] for w in addrs for i in range(len(w))}
    return sorted({v + c for v in internal for c in "01" if v + c not in internal}, key=lambda w: iv(w)[0])

def exps(leaves):
    out = []
    for w in leaves:
        e = len(w) - len(w.rstrip("0")); top = w[:len(w) - e]
        if e > 0 and set(top) <= {"1"}:
            e -= 1
        out.append(e)
    return out

def pos_word(e):
    n = max([k for k in range(1, len(e)) if e[k] > 0], default=0)
    w = "a" * e[0]
    for k in range(1, n + 1):
        w += "b" * e[k]
        if k < n:
            w += "A"
    return w + ("a" * (n - 1) if n >= 1 else "")

def winv(w):
    return "".join(c.swapcase() for c in reversed(w))

def freered(w):
    s = []
    for c in w:
        if s and s[-1] == c.swapcase():
            s.pop()
        else:
            s.append(c)
    return "".join(s)

REF_L, REF_P, REF_M, REF_Q, REF_R = "00", "01", "100", "101", "11"

def pad_last(gap, k):
    gap = list(gap)
    while len(gap) < k:
        w = gap.pop(); gap += [w + "0", w + "1"]
    return gap

def transport(I, J):
    if iv(I)[0] > iv(J)[0]:
        I, J = J, I
    Td = mintree([I, J]); i, j = Td.index(I), Td.index(J)
    gd = [Td[:i], Td[i + 1:j], Td[j + 1:]]
    assert all(len(g) >= 1 for g in gd)                       # (F1): all gaps nonempty
    pr = [pad_last([y], len(x)) for x, y in zip(gd, (REF_L, REF_M, REF_R))]
    Tr = pr[0] + [REF_P] + pr[1] + [REF_Q] + pr[2]
    return I, J, Td, Tr

X0 = [(Fr(0), Fr(0)), (Fr(1, 2), Fr(1, 4)), (Fr(3, 4), Fr(1, 2)), (Fr(1), Fr(1))]
X1 = [(Fr(0), Fr(0)), (Fr(1, 2), Fr(1, 2)), (Fr(3, 4), Fr(5, 8)), (Fr(7, 8), Fr(3, 4)), (Fr(1), Fr(1))]
GEN = {"a": X0, "A": [(y, x) for x, y in X0], "b": X1, "B": [(y, x) for x, y in X1]}

def push(w, c, d):
    """image of [c,d] under the word (rightmost letter first) if affine at every step, else None."""
    for ch in reversed(w):
        ok = False
        for (x0, y0), (x1, y1) in zip(GEN[ch], GEN[ch][1:]):
            if x0 <= c and d <= x1:
                sl = (y1 - y0) / (x1 - x0); c, d = y0 + sl * (c - x0), y0 + sl * (d - x0); ok = True; break
        if not ok:
            return None
    return c, d

def lcp(x, y):
    n = 0
    while n < min(len(x), len(y)) and x[n] == y[n]:
        n += 1
    return n

def analyze(I, J, check=False):
    I, J, Td, Tr = transport(I, J)
    N = len(Td) - 1
    assert N == len(I) + len(J) - lcp(I, J) - 1                              # (F1)
    ed, er = exps(Td), exps(Tr)
    assert sum(ed) == N - len(Td[-1]) and sum(er) == N - len(Tr[-1])         # sum e_k = N - s
    assert ed[-1] == 0 and ed[-2] == 0 and er[-1] == 0 and er[-2] == 0       # last two exponents vanish
    w = freered(pos_word(er) + winv(pos_word(ed)))
    n0 = sum(c in "aA" for c in w); n1 = sum(c in "bB" for c in w)
    e0d, e0r = ed[0], er[0]
    assert e0d == len(I) - len(I.lstrip("0")) and e0r <= 2
    assert n1 <= (N - len(Tr[-1]) - e0r) + (N - len(Td[-1]) - e0d)           # (F2) x1-count
    assert n0 <= e0r + e0d + 2 * N - 6                                       # (F2) x0-count
    g = J.count("0")                                                         # right gap of T_dom
    assert len(Td) - 1 - Td.index(J) == g and len(Tr[-1]) == g + 1             # g_R = #zeros(J), s_ran = g_R+1
    assert all(e == 0 for e in ed[N - g:]) and all(e == 0 for e in er[N - g:]) # exps vanish from J / Q0 on
    assert n0 <= 2 * N - 2 * g + e0d - 2                                     # (F2') x0-count
    assert n1 <= 2 * N - g - 3 - e0d                                         # (F2') x1-count
    if check:
        assert push(w, *iv(I)) == iv(REF_P) and push(w, *iv(J)) == iv(REF_Q)  # it is a transport
    return N, n0, n1, e0d, g

def check_cfp(ntrees=400, seed=3):
    """(R_N -> T) = x0^{e0} x1^{e1} ... x_N^{e_N} (telescoped word pos_word): the word must map
    the k-th leaf of the right vine affinely onto the k-th leaf of T, for every k."""
    import random
    random.seed(seed)
    for _ in range(ntrees):
        n = random.randint(1, 9); T = [""]
        for _ in range(n):
            w = T.pop(random.randrange(len(T))); T += [w + "0", w + "1"]
        T = sorted(T, key=lambda w: iv(w)[0])
        vine = ["1" * i + "0" for i in range(n)] + ["1" * n]
        w = pos_word(exps(T))
        for a, b in zip(vine, T):
            assert push(w, *iv(a)) == iv(b)
    # x_2 = x0^{-1} x1 x0 in the convention hg = h o g maps [3/4,7/8] affinely onto [3/4,13/16]
    assert push("Aba", Fr(3, 4), Fr(7, 8)) == (Fr(3, 4), Fr(13, 16))
    print(f"  CFP normal form (eq. cfp) verified on {ntrees} random trees (<= 9 carets)")

def part_C(mmax=4):
    print("[C] word lengths: exhaustive unit test of (F1),(F2), averages, transports")
    check_cfp()
    for m in range(1, mmax + 1):
        U = ["".join(t) for t in itertools.product("01", repeat=m)]
        tot = [Fr(0)] * 5; W = 0
        terms = []
        for u in U:
            for v, v2 in itertools.permutations(U, 2):
                terms.append((u + "01" + v + "01", u + "01" + v2 + "01", 1))      # siblings
            for w_, w2 in itertools.permutations(U, 2):
                terms.append((w_ + "01", w2 + "01", 1))                           # parents
            for v in U:
                for x in U:
                    if x != u:
                        terms.append((u + "01" + v + "01", x + "01", 2))           # mixed
        cache = {}
        for k, (I, J, wt) in enumerate(terms):
            key = (I, J)
            if key not in cache:
                cache[key] = analyze(I, J, check=(m <= 3 or k % 7 == 0))
            N, n0, n1, e0, g = cache[key]
            tot = [tot[0] + wt * N, tot[1] + wt * n0, tot[2] + wt * n1, tot[3] + wt * e0, tot[4] + wt * g]; W += wt
        Nb, e0b, Mb = closed_form(m)
        avN, av0, av1, ave, avg_ = (t / W for t in tot)
        assert W == 4 * 4 ** m * (2 ** m - 1) and avN == Nb and ave == e0b
        assert av0 <= 2 * Nb + e0b - 4 and av1 <= 2 * Nb - 2
        _, _, gb, Ms = closed_form_sharp(m)
        assert avg_ == gb and max(av0, av1) <= Ms
        print(f"  m={m}: gbar = {avg_} = (3m+4)/4 ok; Mbar' = {float(Ms):.3f} (F2') vs Mbar = {float(Mb):.3f}")
        print(f"  m={m}: weight {W}; avg N = {avN} (closed form ok); avg e0 = {ave} (ok); "
              f"avg n0 = {float(av0):.3f} <= {float(2*Nb+e0b-4):.3f}; avg n1 = {float(av1):.3f} <= {float(2*Nb-2):.3f}")

def part_C_random(seed=7):
    """pointwise (F1),(F2),(F2') and the transport property on random pairs of the layout, large m."""
    import random
    random.seed(seed)
    U = lambda m: "".join(random.choice("01") for _ in range(m))
    for m, ns, nchk in ((8, 3000, 300), (16, 2000, 100), (24, 2000, 50), (40, 1000, 20), (64, 500, 10)):
        for k in range(ns):
            u = U(m); kind = k % 3
            if kind == 0:
                v, v2 = U(m), U(m)
                while v2 == v: v2 = U(m)
                I, J = u + "01" + v + "01", u + "01" + v2 + "01"
            elif kind == 1:
                w2 = U(m)
                while w2 == u: w2 = U(m)
                I, J = u + "01", w2 + "01"
            else:
                x = U(m)
                while x == u: x = U(m)
                I, J = u + "01" + U(m) + "01", x + "01"
            analyze(I, J, check=(k < nchk))
        print(f"  random m={m}: {ns} pairs pass (F1),(F2),(F2'); {nchk} words checked to be transports")

# ---------------------------------------------------------------------------------------------
# [D] main bounds
# ---------------------------------------------------------------------------------------------
def G(kappa, mrange=range(2, 200), sharp=False):
    best = None
    for m in mrange:
        D = 2 ** m
        if kappa <= Fr(4, D):
            continue
        Mb = closed_form_sharp(m)[3] if sharp else closed_form(m)[2]
        c = (kappa - Fr(4, D)) / (4 * (1 - Fr(1, D)) * Mb)
        if best is None or c > best[0]:
            best = (c, m, Mb)
        elif best is not None and m > best[1] + 3:
            break                          # unimodal in m: kappa-4/D saturates, Mbar grows linearly
    return best

def kappa_retraction(k):
    k = Fr(k)
    return (Fr(2) / ((k + 2) * k)) ** 2      # f = -R: displacement 2/(k+2), Lipschitz k

def part_D():
    print("[D] main bounds  h_E(F; x0,x1) >= G(kappa)")
    out = {}
    for name, L in (("unconditional, 1/eps (Thm 1.1)", L_NEW), ("unconditional, 2/eps (Lean as written)", L_OLD)):
        kap = (DELTA / L) ** 2
        c, m, Mb = G(kap)
        out[name] = c
        print(f"  {name}: m={m} (D=2^{m}), Mbar={Mb} = {float(Mb):.4f}, h_E >= {dec(c, 8)}")
        print(f"      exact rational bound has {len(str(c.denominator))}-digit denominator; "
              f"safe stated value: h_E > {dec(c, 3)}")
    for name, L in (("sharpened (F2'), 1/eps", L_NEW),):
        c, m, Mb = G((DELTA / L) ** 2, sharp=True); out[name] = c
        print(f"  {name}: m={m}, Mbar'={float(Mb):.4f}, h_E >= {dec(c, 8)}  (gain {dec(c / out['unconditional, 1/eps (Thm 1.1)'], 5)})")
    r = out["unconditional, 1/eps (Thm 1.1)"] / out["unconditional, 2/eps (Lean as written)"]
    print(f"  ratio new/old = {dec(r, 6)}")
    print("  conditional g(k) (k-Lipschitz retraction B -> S of the Hilbert ball):")
    table = {}
    for k in (Fr(3), Fr(9, 2), Fr(458, 100), 5, 6, 8, 10, 15, 20, 25, Fr(2899, 100), 29, Fr(3226, 100), 33, 40, 50, 100):
        c, m, Mb = G(kappa_retraction(k))
        table[k] = (c, m, Mb)
        cs, ms, Ms = G(kappa_retraction(k), sharp=True)
        table[k] = (c, m, Mb, cs)
        print(f"    k={float(k):7.2f}: kappa={dec(kappa_retraction(k), 6):>12s}  m={m:3d}  Mbar={float(Mb):8.3f}  g(k) = {dec(c, 6)}   [sharp: m={ms}, {dec(cs, 6)}]")
    return out, table

# ---------------------------------------------------------------------------------------------
# [E] ceilings, ledger, spectral
# ---------------------------------------------------------------------------------------------
def part_E(out):
    print("[E] framework ceilings / loss ledger / spectral corollary")
    c, m, Mb = G(Fr(1, 16))
    print(f"  kappa <= 1/16 for every Lipschitz f:B->B (Goebel: delta <= 1-1/L).  Reference-pair ceiling G(1/16) = "
          f"{dec(c, 6)} at m={m}; S-form ceiling kappa/4 <= 1/64")
    # ledger for the unconditional bound (log2 of each factor)
    h = out["unconditional, 1/eps (Thm 1.1)"]
    kap = (DELTA / L_NEW) ** 2
    cc, mm, MM = G(kap)
    lg = lambda x: math.log2(float(x))
    print("  ledger (log2):")
    print(f"    delta^2 = 1/4                        : {lg(DELTA**2):9.3f}")
    print(f"    1/eps^2 = (4/rho)^2                  : {-2*lg(1/EPS):9.3f}   (rho = 1/2048)")
    print(f"    1/Lip(G)^2                           : {-2*lg(LIPG):9.3f}   (Lip(G) ~ 7/(4 rho))")
    print(f"    kappa                                : {lg(kap):9.3f}")
    print(f"    (1-4/(D kappa))/(4(1-1/D))           : {lg((kap-Fr(4,2**mm))/kap/(4*(1-Fr(1,2**mm)))):9.3f}")
    print(f"    1/Mbar                               : {-lg(MM):9.3f}   (Mbar = {float(MM):.3f}, m = {mm})")
    print(f"    h_E                                  : {lg(h):9.3f}")
    # spectral: Mohar 1988, lambda_0(Delta) >= d - sqrt(d^2 - h^2) >= h^2/(2d), d = 4
    for name, hv in out.items():
        gap = hv ** 2 / 8
        print(f"  {name}: Laplacian bottom of spectrum >= h^2/8 = {dec(gap, 4)}; SRW spectral radius <= 1 - h^2/32 = 1 - {dec(hv**2/32, 4)}")

if __name__ == "__main__":
    fast = "--fast" in sys.argv
    if not fast:
        part_A()
    part_B()
    part_C(3 if fast else 4)
    if not fast:
        part_C_random()
    out, table = part_D()
    part_E(out)
    print("certificate complete")
