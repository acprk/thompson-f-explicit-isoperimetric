#!/usr/bin/env python3
"""
Explicit constants (with the Banach-space 2/eps normalization; cert.py uses 1/eps), exact rationals + arb for the final decimals.

Notation: kappa = (delta/L)^2 for the displacement map f:B->B (||f(x)-x||>=delta, Lip f <= L).
  Paper Prop avg:boundary:      max_{h in S} |hA(+)A|/|A| >= (kappa - 4/D)/(4(1-1/D)).
  Reference-pair route ({x0,x1}): |x0A(+)A|/|A| + |x1A(+)A|/|A| >= (kappa - 4/D)/(4(1-1/D) Mbar(m)),
      D = 2^m, Mbar(m) = max(2 Nbar + e0bar - 4, 2 Nbar - 2) from s4_wordlength_R.closed_form.
  Max-version (paper's max over S, no averaging): Mmax(m) = 7r - 7, r = m+2.
Inputs for f:
  (P) paper/Lean map made explicit by s1: delta = 1/2, L = 59242240.
  (K) -R for a 29-Lipschitz retraction R: B -> S of the Hilbert ball (conditional: k0(H) <= 28.99 is attributed
      to Baronti-Casini-Franchetti 2003 via Mahesh arXiv:2609.13342; not checked against the original): delta = 2/(2+29) = 2/31, L = 29.
"""
import sys
from fractions import Fraction as Fr
sys.path.insert(0, ".")
from s4_wordlength_R import closed_form
from flint import arb, ctx
ctx.prec = 120

def paper_S_level(kappa):
    # sup over integer D of (kappa-4/D)/(4(1-1/D)); increasing in D -> limit kappa/4
    return kappa / 4

def Mbar(m):
    Nb, e0 = closed_form(m)
    return max(2 * Nb + e0 - 4, 2 * Nb - 2)

def best_R(kappa, mrange, avg=True):
    best = None
    for m in mrange:
        D = 2 ** m
        num = kappa - Fr(4, D)
        if num <= 0:
            continue
        M = Mbar(m) if avg else Fr(7 * (m + 2) - 7)
        c = num / (4 * (1 - Fr(1, D)) * M)
        if best is None or c > best[0]:
            best = (c, m, M)
    return best

def show(name, x):
    print(f"  {name:58s} = {arb(x.numerator) / arb(x.denominator)}")

if __name__ == "__main__":
    LP = Fr(59242240); kP = Fr(1, 4) / LP ** 2
    kK = (Fr(2, 31) / 29) ** 2
    print("kappa_paper =", float(kP), "  kappa_K =", float(kK), "  ratio", float(kK / kP))
    print("[S-level, paper's own S (transports), D -> infinity]")
    show("paper explicit  kappa_P/4", paper_S_level(kP))
    show("paper framework with f=-R  kappa_K/4", paper_S_level(kK))
    print("[{x0,x1} edge-isoperimetric constant  h_E = inf (|x0A(+)A|+|x1A(+)A|)/|A|]")
    for label, k, mr in (("paper map, max-version", kP, range(40, 80)),
                         ("paper map, averaged (reference-pair route)", kP, range(40, 80)),
                         ("f=-R, max-version", kK, range(10, 40)),
                         ("f=-R, averaged (reference-pair route)  [MAIN]", kK, range(10, 40))):
        c, m, M = best_R(k, mr, avg=("averaged" in label))
        print(f" {label}: m={m} (D=2^{m}, r={m+2}), Mbar={float(M):.4f}")
        show("lower bound", c)
    # ceilings
    print("[framework ceilings]")
    # (i) any f: delta <= 1 - 1/L  =>  delta/L <= max (1-1/L)/L = 1/4
    show("kappa ceiling (any f, self-proved delta<=1-1/L)", Fr(1, 16))
    # Hilbert bound psi_H(k) <= sqrt(1-1/k) k/(k+1)  (Goebel/Bolibok, cited):  kappa <= max (1-1/k)/(k+1)^2
    best = max(((1 - 1 / k) / (k + 1) ** 2, k) for k in [1 + i / 10000 for i in range(1, 100000)])
    print(f"  kappa ceiling from psi_H(k)<=sqrt(1-1/k)k/(k+1): {best[0]:.6f} at k={best[1]:.4f}")
    # (ii) reference-pair word ceiling: slope argument max(n0,n1) >= (n0+n1)/2 >= (d_I - 2)/2 averaged
    for m in (21, 22, 23):
        r = m + 2
        # weighted average of the deeper interval depth over S(2r),P(r),M(2r): (2r + r + 2*2r)/4
        lowM = Fr(2 * r + r + 2 * 2 * r, 4) / 2 - 1
        print(f"  m={m}: reference-pair slope floor  max(n0,n1)avg >= {float(lowM):.2f}  vs our Mbar {float(Mbar(m)):.2f}")
