#!/usr/bin/env python3
"""Independent re-implementation (independent of cert.py / v2_words.py data structures) of (F1), (F2'),
the closed-form averages and the transport property.

Elements of F are exact PL maps on dyadics, stored as integer pairs over a common power of two.
Trees are built as explicit node objects; exponents are computed by walking parent pointers.
The word for (T_dom -> T_ran) is P(T_ran) P(T_dom)^{-1}, freely reduced; we check that the *whole
PL map* of the word equals the tree-pair map leaf by leaf (stronger than checking I and J only).
"""
import sys, random, itertools
from fractions import Fraction as Fr
from multiprocessing import Pool

# ---------- PL maps: list of (x,y) Fractions, increasing ----------
def pl_eval(f, x):
    lo, hi = 0, len(f) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if f[mid][0] <= x: lo = mid
        else: hi = mid
    (x0, y0), (x1, y1) = f[lo], f[hi]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)

def pl_is_affine_on(f, a, b):
    return not any(a < x < b for x, _ in f)

X0 = [(Fr(0), Fr(0)), (Fr(1, 2), Fr(1, 4)), (Fr(3, 4), Fr(1, 2)), (Fr(1), Fr(1))]
X1 = [(Fr(0), Fr(0)), (Fr(1, 2), Fr(1, 2)), (Fr(3, 4), Fr(5, 8)), (Fr(7, 8), Fr(3, 4)), (Fr(1), Fr(1))]
inv = lambda f: [(y, x) for x, y in f]
GEN = {'x': X0, 'X': inv(X0), 'y': X1, 'Y': inv(X1)}

def word_maps_interval(word, a, b):
    """apply word (rightmost letter first); return image of [a,b] if affine throughout, else None"""
    for c in reversed(word):
        f = GEN[c]
        if not pl_is_affine_on(f, a, b):
            return None
        a, b = pl_eval(f, a), pl_eval(f, b)
    return a, b

def interval(addr):
    k = int(addr, 2) if addr else 0
    return Fr(k, 1 << len(addr)), Fr(k + 1, 1 << len(addr))

# ---------- trees as explicit nodes ----------
class Node:
    __slots__ = ("addr", "parent", "kids")
    def __init__(s, addr, parent):
        s.addr, s.parent, s.kids = addr, parent, None

def build_tree(leafset):
    """leafset: set of binary words forming a complete prefix code. Returns ordered leaf Nodes."""
    root = Node("", None)
    leaves = []
    def rec(nd):
        if nd.addr in leafset:
            leaves.append(nd); return
        nd.kids = (Node(nd.addr + "0", nd), Node(nd.addr + "1", nd))
        rec(nd.kids[0]); rec(nd.kids[1])
    rec(root)
    return leaves

def on_right_spine(nd):
    return all(ch == "1" for ch in nd.addr)

def exponent(leaf):
    """CFP: length of maximal chain of left edges ending at leaf, not counting edges hanging from
    a right-spine vertex."""
    e, nd = 0, leaf
    while nd.parent is not None and nd.addr[-1] == "0":
        if on_right_spine(nd.parent):
            break
        e += 1; nd = nd.parent
    return e

def nf_word(exps):
    """x_0^{e0} x_1^{e1} ... x_N^{eN}, x_k = x0^{-(k-1)} x1 x0^{k-1}; build literally then reduce."""
    w = "x" * exps[0]
    for k in range(1, len(exps)):
        if exps[k]:
            w += "X" * (k - 1) + "y" * exps[k] + "x" * (k - 1)
    return reduce_word(w)

def reduce_word(w):
    st = []
    for c in w:
        if st and st[-1] == c.swapcase(): st.pop()
        else: st.append(c)
    return "".join(st)

def winv(w): return "".join(c.swapcase() for c in reversed(w))

# ---------- the transport ----------
def leaves_of_min_tree(I, J):
    pref = {w[:i] for w in (I, J) for i in range(len(w))}
    return {p + c for p in pref for c in "01" if p + c not in pref}

def split_last(gap, k):
    gap = list(gap)
    while len(gap) < k:
        last = gap.pop(); gap.extend([last + "0", last + "1"])
    return gap

def analyse(I, J, full_check):
    if interval(I)[0] > interval(J)[0]:
        I, J = J, I
    Ld = build_tree(leaves_of_min_tree(I, J))
    addrs = [nd.addr for nd in Ld]
    i, j = addrs.index(I), addrs.index(J)
    gaps = (addrs[:i], addrs[i + 1:j], addrs[j + 1:])
    assert min(map(len, gaps)) >= 1
    ran = split_last(["00"], len(gaps[0])) + ["01"] + split_last(["100"], len(gaps[1])) + ["101"] + \
          split_last(["11"], len(gaps[2]))
    Lr = build_tree(set(ran))
    assert [nd.addr for nd in Lr] == ran
    N = len(Ld) - 1
    lcp = next((k for k in range(min(len(I), len(J))) if I[k] != J[k]), min(len(I), len(J)))
    assert N == len(I) + len(J) - lcp - 1                                    # (F1)
    ed = [exponent(nd) for nd in Ld]; er = [exponent(nd) for nd in Lr]
    w = reduce_word(nf_word(er) + winv(nf_word(ed)))
    n0 = sum(c in "xX" for c in w); n1 = sum(c in "yY" for c in w)
    g = J.count("0"); e0d = ed[0]
    ok_sharp = (n0 <= 2 * N - 2 * g + e0d - 2) and (n1 <= 2 * N - g - 3 - e0d)
    ok_old = (n0 <= ed[0] + er[0] + 2 * N - 6) and \
             (n1 <= (N - len(ran[-1]) - er[0]) + (N - len(addrs[-1]) - ed[0]))
    # slack of the sharp bounds (0 means the bound is attained)
    slack = (2 * N - 2 * g + e0d - 2 - n0, 2 * N - g - 3 - e0d - n1)
    if full_check:
        # whole map check: every leaf of T_dom goes affinely onto the matching leaf of T_ran
        for a, b in zip(addrs, ran):
            assert word_maps_interval(w, *interval(a)) == interval(b), (I, J, a, b)
    return N, n0, n1, e0d, g, ok_sharp, ok_old, slack

def layout_pairs(m):
    U = ["".join(t) for t in itertools.product("01", repeat=m)]
    for u in U:
        for v, v2 in itertools.permutations(U, 2):
            yield (u + "01" + v + "01", u + "01" + v2 + "01", 1)
        for w1, w2 in itertools.permutations(U, 2):
            yield (w1 + "01", w2 + "01", 1)
        for v in U:
            for x in U:
                if x != u:
                    yield (u + "01" + v + "01", x + "01", 2)

def closed(m):
    D = 2 ** m; r = m + 2
    lam = Fr(D - 1 - m, D - 1)
    Nb = Fr(11 * r, 4) - 1 - lam
    a = lambda k: D // 2 ** k if k <= m else (1 if k == m + 1 else 0)
    EP = sum(1 - Fr((D - a(k)) * (D - a(k) - 1), D * (D - 1)) for k in range(1, m + 2))
    e0 = (1 + 3 * EP) / 4
    gb = Fr(3 * m + 4, 4)
    M = max(2 * Nb + e0 - 4, 2 * Nb - 2)
    Ms = max(2 * Nb - 2 * gb + e0 - 2, 2 * Nb - gb - 3 - e0)
    return Nb, e0, gb, M, Ms

def work(args):
    I, J, full = args
    return analyse(I, J, full)

def exhaustive(m, pool):
    seen = {}
    terms = list(layout_pairs(m))
    keys = {}
    for I, J, wt in terms:
        k = (I, J) if interval(I)[0] < interval(J)[0] else (J, I)
        keys[k] = None
    klist = list(keys)
    full = m <= 4
    res = pool.map(work, [(I, J, full or (idx % 50 == 0)) for idx, (I, J) in enumerate(klist)], chunksize=500)
    R = dict(zip(klist, res))
    tot = [Fr(0)] * 5; W = 0; bad = 0; minslack = [10 ** 9, 10 ** 9]
    for I, J, wt in terms:
        k = (I, J) if interval(I)[0] < interval(J)[0] else (J, I)
        N, n0, n1, e0, g, oks, oko, sl = R[k]
        if not (oks and oko): bad += 1
        minslack = [min(minslack[0], sl[0]), min(minslack[1], sl[1])]
        for t, val in enumerate((N, n0, n1, e0, g)):
            tot[t] += wt * val
        W += wt
    av = [t / W for t in tot]
    Nb, e0b, gb, M, Ms = closed(m)
    print(f"m={m}: {len(klist)} distinct pairs, weight {W}, failures {bad}, min slack (n0,n1) {minslack}")
    print(f"   avg N={av[0]} closed {Nb} {'OK' if av[0]==Nb else 'MISMATCH'}; avg e0={av[3]} closed {e0b} "
          f"{'OK' if av[3]==e0b else 'MISMATCH'}; avg g={av[4]} closed {gb} {'OK' if av[4]==gb else 'MISMATCH'}")
    print(f"   exact avg n0={float(av[1]):.4f} n1={float(av[2]):.4f}; Mbar={float(M):.4f} Mbar'={float(Ms):.4f}"
          f"  admissible: {max(av[1],av[2]) <= Ms}")

def rnd(m, cnt, rng):
    U = lambda: "".join(rng.choice("01") for _ in range(m))
    out = []
    for k in range(cnt):
        u = U(); kind = k % 3
        if kind == 0:
            v = U(); v2 = U()
            while v2 == v: v2 = U()
            out.append((u + "01" + v + "01", u + "01" + v2 + "01"))
        elif kind == 1:
            w = U()
            while w == u: w = U()
            out.append((u + "01", w + "01"))
        else:
            x = U()
            while x == u: x = U()
            out.append((u + "01" + U() + "01", x + "01"))
    return out

if __name__ == "__main__":
    mmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    with Pool(4) as pool:
        for m in range(1, mmax + 1):
            exhaustive(m, pool)
        rng = random.Random(20261009)
        for m, cnt in ((7, 20000), (10, 20000), (20, 10000), (40, 5000), (59, 5000), (64, 3000), (100, 1000)):
            ps = rnd(m, cnt, rng)
            res = pool.map(work, [(I, J, idx < 200) for idx, (I, J) in enumerate(ps)], chunksize=100)
            bad = sum(1 for r in res if not (r[5] and r[6]))
            ms = (min(r[7][0] for r in res), min(r[7][1] for r in res))
            print(f"random m={m}: {cnt} pairs, failures {bad}, min slack {ms}, 200 full-map checks passed")
    for m in (24, 59):
        Nb, e0b, gb, M, Ms = closed(m)
        print(f"closed m={m}: Nbar={float(Nb):.6f} e0bar={float(e0b):.6f} gbar={gb} Mbar={float(M):.6f} Mbar'={float(Ms):.6f}")
