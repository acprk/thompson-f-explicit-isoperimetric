#!/usr/bin/env python3
"""
Rigorous word-length bookkeeping for the reference-pair route.

Layout (all addresses over {0,1}, '0' = left half):
  parents     P_u = u+'01'          u in {0,1}^m,   D = 2^m, depth r = m+2
  descendants P_u.P_v = u01v01      (I.J = s_I(J) is address concatenation)
  reference   P0 = '01' = [1/4,1/2],  Q0 = '101' = [5/8,3/4]
For each separated pair I<J used in the variance expansion we transport (I,J) onto (P0,Q0) by
  phi_{I,J} = element(T_dom, T_ran),  T_dom = minimal tree with leaves I,J,
  T_ran = {00,01,100,101,11} padded gap-by-gap ("split the last leaf of the gap") to T_dom's gap sizes.
Facts proved in the note and CHECKED here exhaustively for small m:
  (F1) every gap of T_dom is nonempty (I internal, separated) and the reference tree has gaps (1,1,1),
       so no padding on the domain side and N(phi) = d_I + d_J - lcp(I,J) - 1;
  (F2) the explicit (freely reduced) normal-form word w(phi) = P(T_ran) P(T_dom)^{-1} (CFP normal form, telescoped
       through x_k = x0^{-(k-1)} x1 x0^{k-1}) satisfies
         n1(w) <= (N - s_ran - e0_ran) + (N - s_dom - e0_dom),
         n0(w) <= e0_ran + e0_dom + 2N - 6,     e0_ran <= 2,  e0_dom = #leading zeros of I,
       where n0/n1 = number of letters x0^{+-1} / x1^{+-1}.
Then |phi A (+) A| <= n0 |x0 A (+) A| + n1 |x1 A (+) A|.
"""
import sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, ".")
from thompson import (minimal_tree, leaves_sorted, _gaps, _pad, element, nf_word, eval_word,
                      letter_counts, interval, check_transport, free_reduce)

REF_TREE = {"00", "01", "100", "101", "11"}
P0, Q0 = "01", "101"

def lcp(a, b):
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    return n

def lz(a):
    n = 0
    while n < len(a) and a[n] == "0":
        n += 1
    return n

def phi_trees(I, J):
    Td = minimal_tree([I, J])
    gd = _gaps(Td, I, J)
    gr = _gaps(REF_TREE, P0, Q0)
    assert all(len(x) >= 1 for x in gd) and all(len(y) == 1 for y in gr)
    newr = [_pad(y, len(x), "last") for x, y in zip(gd, gr)]
    Tr = set(newr[0] + [P0] + newr[1] + [Q0] + newr[2])
    return Td, Tr

def depth_leftmost(T):
    return len(leaves_sorted(T)[0])

def depth_rightmost(T):
    return len(leaves_sorted(T)[-1])

def phi_data(I, J, verify=False):
    if interval(I)[0] > interval(J)[0]:
        I, J = J, I
    Td, Tr = phi_trees(I, J)
    N = len(Td) - 1
    assert N == len(I) + len(J) - lcp(I, J) - 1                 # (F1)
    w = word_unreduced(Td, Tr)
    n0, n1 = letter_counts(w)
    e0d, e0r = depth_leftmost(Td) - 1, depth_leftmost(Tr) - 1
    sd, sr = depth_rightmost(Td), depth_rightmost(Tr)
    assert e0d == lz(I) and e0r <= 2                               # (F2) side facts
    assert n1 <= (N - sr - e0r) + (N - sd - e0d)                   # (F2)
    assert n0 <= e0r + e0d + 2 * N - 6                             # (F2)
    if verify:
        assert check_transport(w, I, J, P0, Q0)
    return N, n0, n1, e0d

def word_unreduced(Tdom, Tran):
    from thompson import positive_word_from_tree
    pos, neg = positive_word_from_tree(Tran), positive_word_from_tree(Tdom)
    def pos_part(p):
        n = max([k for k, e in enumerate(p) if e > 0 and k >= 1], default=0)
        w = "a" * p[0]
        for k in range(1, n + 1):
            w += "b" * p[k]
            if k < n:
                w += "A"
        if n >= 1:
            w += "a" * (n - 1)
        return w
    inv = lambda w: "".join(c.swapcase() for c in reversed(w))
    # free reduction (same group element); only cancels letters, so counts below are upper-bounded
    return free_reduce(pos_part(pos) + inv(pos_part(neg)))

def exhaustive(m, verify_frac=1.0):
    U = ["".join(t) for t in itertools.product("01", repeat=m)]
    D = len(U)
    tot = {"n0": Fr(0), "n1": Fr(0), "N": Fr(0), "e0": Fr(0), "w": Fr(0)}
    cnt = 0
    k = 0
    for u in U:
        terms = []
        for v, v2 in itertools.permutations(U, 2):           # (+S) ordered sibling pairs, weight 1
            terms.append((u + "01" + v + "01", u + "01" + v2 + "01", 1))
        for w_, w2 in itertools.permutations(U, 2):          # (+P) ordered parent pairs, weight 1
            terms.append((w_ + "01", w2 + "01", 1))
        for v in U:                                            # (-M) mixed, weight 2
            for x in U:
                if x != u:
                    terms.append((u + "01" + v + "01", x + "01", 2))
        for I, J, wt in terms:
            k += 1
            N, n0, n1, e0d = phi_data(I, J, verify=(k % int(1 / verify_frac) == 0))
            tot["n0"] += wt * n0; tot["n1"] += wt * n1; tot["N"] += wt * N; tot["e0"] += wt * e0d
            cnt += wt
    return {key: val / cnt for key, val in tot.items()}, cnt, 4 * D * D * (D - 1)

def closed_form(m):
    """exact weighted averages (weights: S 1, P 1, M 2 over all rows u) of N and e0_dom."""
    D = 2 ** m; r = m + 2
    lam = Fr(D - m - 1, D - 1)              # E[lcp(v,v') | v != v']  (same for lcp(u,x), x != u)
    Nbar = (Fr(3 * r - 1) - lam + Fr(2 * r - 1) - lam + 2 * (Fr(3 * r - 1) - lam)) / 4
    # e0_dom = leading zeros of the lower interval
    def lz_ge(k):  # #{u : lz(u01) >= k}
        if k <= m: return D // 2 ** k
        if k == m + 1: return 1
        return 0
    E_lz = sum(Fr(lz_ge(k), D) for k in range(1, m + 2))                       # = 1
    # (+P): lower = min(w,w'); lz(min) = max(lz(w),lz(w')); pairs w != w'
    def max_ge(k):
        a = lz_ge(k); nb = D - a
        return Fr(D * (D - 1) - nb * (nb - 1), D * (D - 1))
    E_P = sum(max_ge(k) for k in range(1, m + 2))
    # (-M): lower is x01 if x<u else u01..; average over u, x != u of lz(lower)
    #       x<u  <=> lower = x; lz(min(u,x)) = max(lz(u),lz(x)) as well
    E_M = E_P
    E_S = E_lz
    e0bar = (E_S + E_P + 2 * E_M) / 4
    return Nbar, e0bar

if __name__ == "__main__":
    for m in (1, 2, 3, 4):
        avg, cnt, expect = exhaustive(m, verify_frac=1.0 if m <= 3 else 0.2)
        Nbar, e0bar = closed_form(m)
        print(f"m={m}: weight={cnt} (=4D^2(D-1): {cnt == expect})  exact avg N={avg['N']} closed={Nbar}  "
              f"avg e0={avg['e0']} closed={e0bar}  avg n0={float(avg['n0']):.3f} avg n1={float(avg['n1']):.3f}  "
              f"bound 2N+e0+2-6={float(2 * Nbar + e0bar - 4):.3f}")
        assert avg["N"] == Nbar and avg["e0"] == e0bar
        assert avg["n0"] <= 2 * Nbar + e0bar - 4 and avg["n1"] <= 2 * Nbar + e0bar - 4
    print("all pointwise facts (F1),(F2) and transport properties verified; closed forms exact.")
