#!/usr/bin/env python3
"""
Exact Thompson's group F toolkit (Fractions only, no floating point).

Elements = increasing PL homeomorphisms of [0,1] stored as breakpoint lists (xs, ys).
Multiplication convention of the paper: (h g)(t) = h(g(t)).
Generators (Cannon-Floyd-Parry): x0 = [0,1/2]->[0,1/4], [1/2,3/4]->[1/4,1/2], [3/4,1]->[1/2,1];
x1 = id on [0,1/2] and a half-scale copy of x0 on [1/2,1].  With h∘g composition,
x_{n+1} = x0^{-1} x_n x0, i.e. x_n = x0^{-(n-1)} x1 x0^{n-1}  (checked in selftest()).

Tree-pair diagrams: a finite rooted binary tree is given by its set of leaf addresses
(strings over {'0','1'}, '0' = left half); leaf w <-> dyadic interval [0.w, 0.w + 2^-|w|].
The element (Tdom -> Tran) maps the k-th leaf interval of Tdom affinely onto the k-th of Tran.
"""
from fractions import Fraction as Fr
from functools import reduce

# ---------- PL maps ----------
class PL:
    __slots__ = ("xs", "ys")
    def __init__(self, xs, ys):
        self.xs = list(xs); self.ys = list(ys)
        self._normalize()
    def _normalize(self):
        xs, ys = self.xs, self.ys
        keep_x, keep_y = [xs[0]], [ys[0]]
        for i in range(1, len(xs) - 1):
            s1 = (ys[i] - keep_y[-1]) / (xs[i] - keep_x[-1])
            s2 = (ys[i + 1] - ys[i]) / (xs[i + 1] - xs[i])
            if s1 != s2:
                keep_x.append(xs[i]); keep_y.append(ys[i])
        keep_x.append(xs[-1]); keep_y.append(ys[-1])
        self.xs, self.ys = keep_x, keep_y
    def __call__(self, t):
        xs, ys = self.xs, self.ys
        lo, hi = 0, len(xs) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if xs[mid] <= t: lo = mid
            else: hi = mid
        return ys[lo] + (ys[lo + 1] - ys[lo]) * (t - xs[lo]) / (xs[lo + 1] - xs[lo])
    def inv(self):
        return PL(self.ys, self.xs)
    def __matmul__(self, g):  # self ∘ g
        pts = set(g.xs) | {g.inv()(x) for x in self.xs}
        xs = sorted(pts)
        return PL(xs, [self(g(x)) for x in xs])
    def __eq__(self, o):
        return self.xs == o.xs and self.ys == o.ys
    def is_identity(self):
        return self.xs == [Fr(0), Fr(1)]
    def affine_on(self, a, b, c, d):
        """True iff self maps [a,b] affinely onto [c,d]."""
        if self(a) != c or self(b) != d:
            return False
        return not any(a < x < b for x in self.xs)
    def __repr__(self):
        return "PL(" + ", ".join(f"{x}->{y}" for x, y in zip(self.xs, self.ys)) + ")"

ID = PL([Fr(0), Fr(1)], [Fr(0), Fr(1)])
X0 = PL([Fr(0), Fr(1, 2), Fr(3, 4), Fr(1)], [Fr(0), Fr(1, 4), Fr(1, 2), Fr(1)])
X1 = PL([Fr(0), Fr(1, 2), Fr(3, 4), Fr(7, 8), Fr(1)], [Fr(0), Fr(1, 2), Fr(5, 8), Fr(3, 4), Fr(1)])
GEN = {"a": X0, "A": X0.inv(), "b": X1, "B": X1.inv()}  # a=x0, A=x0^-1, b=x1, B=x1^-1

def free_reduce(w):
    out = []
    for c in w:
        if out and out[-1] == c.swapcase():
            out.pop()
        else:
            out.append(c)
    return "".join(out)

def eval_word(w):
    """word s1 s2 ... sk  ->  s1 ∘ s2 ∘ ... ∘ sk"""
    g = ID
    for c in reversed(w):
        g = GEN[c] @ g
    return g

def letter_counts(w):
    return sum(c in "aA" for c in w), sum(c in "bB" for c in w)

# ---------- trees ----------
def interval(addr):
    n = len(addr)
    k = int(addr, 2) if n else 0
    return Fr(k, 2 ** n), Fr(k + 1, 2 ** n)

def leaves_sorted(T):
    return sorted(T, key=lambda w: interval(w)[0])

def is_tree(T):
    """leaf set of a full binary tree: dyadic intervals tiling [0,1]."""
    L = leaves_sorted(T)
    if interval(L[0])[0] != 0 or interval(L[-1])[1] != 1:
        return False
    return all(interval(L[i])[1] == interval(L[i + 1])[0] for i in range(len(L) - 1))

def minimal_tree(addrs):
    """smallest tree having every address in addrs as a leaf (addrs pairwise non-nested)."""
    internal = {""}
    for w in addrs:
        for i in range(len(w)):
            internal.add(w[:i])
    leaves = set()
    for v in internal:
        for c in "01":
            if v + c not in internal:
                leaves.add(v + c)
    if addrs == [""] or (len(addrs) == 1 and addrs[0] == ""):
        return {""}
    return leaves

def element(Tdom, Tran):
    A, B = leaves_sorted(Tdom), leaves_sorted(Tran)
    assert len(A) == len(B)
    xs = [interval(w)[0] for w in A] + [Fr(1)]
    ys = [interval(w)[0] for w in B] + [Fr(1)]
    return PL(xs, ys)

def right_vine(n):
    """right vine with n carets: leaves 0, 10, 110, ..., 1^{n-1}0, 1^n."""
    return {"1" * i + "0" for i in range(n)} | {"1" * n}

def leaf_exponents(T):
    """CFP leaf exponent: for leaf k, the length of the maximal path of left edges ending at
    leaf k whose top vertex is NOT on the right spine (the right side of the tree)."""
    L = leaves_sorted(T)
    exps = []
    for w in L:
        e = 0; v = w
        while v and v[-1] == "0":
            v = v[:-1]; e += 1
        # v is the top vertex of the maximal left path; it lies on the right spine iff v is all '1'
        if set(v) <= {"1"}:
            e = 0  # path reaches the right side -> exponent 0 ... (except leftmost leaf, see below)
        exps.append(e)
    return L, exps

def positive_word_from_tree(T):
    """Candidate positive normal form for the element (T -> right vine) or (right vine -> T);
    the convention is fixed by selftest()."""
    L = leaves_sorted(T)
    exps = []
    for w in L:
        e = 0; v = w
        while v and v[-1] == "0":
            v = v[:-1]; e += 1
        # CFP: maximal left path from the leaf that does NOT reach the right side:
        # if its top vertex v lies on the right spine (v in 1*), drop the last edge.
        if e > 0 and set(v) <= {"1"}:
            e -= 1
        exps.append(e)
    return exps

def nf_word(Tdom, Tran):
    """word in {x0,x1} for element(Tdom, Tran) = P(Tran) P(Tdom)^{-1}, P(T)=element(vine,T)."""
    return word_from_exponents(positive_word_from_tree(Tran), positive_word_from_tree(Tdom))

def word_from_exponents(pos, neg):
    """x_0^{p0} x_1^{p1}...x_n^{pn} x_m^{-q_m} ... x_0^{-q0} written in {x0,x1} using
    x_k = x0^{-(k-1)} x1 x0^{k-1}, telescoped."""
    def pos_part(p):
        n = max([k for k, e in enumerate(p) if e > 0 and k >= 1], default=0)
        w = "a" * p[0]
        for k in range(1, n + 1):
            w += "b" * p[k]
            if k < n: w += "A"
        w += "a" * (n - 1) if n >= 1 else ""
        return w
    def inverse(w):
        return "".join(c.swapcase() for c in reversed(w))
    return free_reduce(pos_part(pos) + inverse(pos_part(neg)))

def selftest():
    # x_{n+1} = x0^{-1} x_n x0
    x2 = eval_word("Aba")
    assert x2.xs[:2] == [Fr(0), Fr(3, 4)] and x2(Fr(7, 8)) == Fr(13, 16), x2
    # tree pair of x0 : ({0, 10, 11} -> {00, 01, 1})
    assert element({"0", "10", "11"}, {"00", "01", "1"}) == X0
    assert element({"0", "10", "110", "111"}, {"0", "100", "101", "11"}) == X1
    return True

if __name__ == "__main__":
    print("selftest", selftest())

# ---------- transports (Lemma dyad:transport, constructive) ----------
def _gaps(T, c, d):
    L = leaves_sorted(T)
    i, j = L.index(c), L.index(d)
    return L[:i], L[i + 1:j], L[j + 1:]

def _pad(gap, k, rule):
    gap = list(gap)
    while len(gap) < k:
        if rule == "shallow":
            idx = min(range(len(gap)), key=lambda t: (len(gap[t]), -t))
        elif rule == "last":
            idx = len(gap) - 1
        elif rule == "first":
            idx = 0
        else:
            raise ValueError(rule)
        w = gap.pop(idx)
        gap[idx:idx] = [w + "0", w + "1"]
    return gap

def transport_pair(c, d, a, b, rule="shallow"):
    """tree pair (Tdom, Tran) of an element carrying interval(c) affinely onto interval(a)
    and interval(d) onto interval(b); requires c<d, a<b separated internal basic intervals."""
    Td, Tr = minimal_tree([c, d]), minimal_tree([a, b])
    gd, gr = _gaps(Td, c, d), _gaps(Tr, a, b)
    assert all((len(x) > 0) == (len(y) > 0) for x, y in zip(gd, gr)), "inconsistent empty gaps"
    newd, newr = [], []
    for x, y in zip(gd, gr):
        k = max(len(x), len(y))
        newd.append(_pad(x, k, rule)); newr.append(_pad(y, k, rule))
    Dl = newd[0] + [c] + newd[1] + [d] + newd[2]
    Rl = newr[0] + [a] + newr[1] + [b] + newr[2]
    return set(Dl), set(Rl)

def check_transport(w, c, d, a, b):
    g = eval_word(w)
    return g.affine_on(*interval(c), *interval(a)) and g.affine_on(*interval(d), *interval(b))

# ---------- suffix-sharing transports ----------
def lo_end(w):
    return interval(w)[0]

def disjoint(p, q):
    return not (p.startswith(q) or q.startswith(p))

def gap_signature(Tleaves, t1, t2):
    g = _gaps(Tleaves, t1, t2)
    return tuple(len(x) > 0 for x in g)

def suffix_transport(c, d, a, b, maxcut=None):
    """Search leaves tc|c, td|d, ta|a, tb|b (prefixes) with c=tc+s1, a=ta+s1, d=td+s2, b=tb+s2,
    giving a valid tree pair (consistent gap emptiness); returns the one with the fewest carets.
    c<d and a<b required (interval order)."""
    best = None
    # common suffix lengths
    def cs(x, y):
        n = 0
        while n < min(len(x), len(y)) and x[-1 - n] == y[-1 - n]:
            n += 1
        return n
    s1max, s2max = cs(c, a), cs(d, b)
    for s1 in range(0, s1max + 1):
        tc, ta = c[:len(c) - s1], a[:len(a) - s1]
        if tc == "" or ta == "":
            continue
        for s2 in range(0, s2max + 1):
            td, tb = d[:len(d) - s2], b[:len(b) - s2]
            if td == "" or tb == "":
                continue
            if not (disjoint(tc, td) and disjoint(ta, tb)):
                continue
            if not (lo_end(tc) < lo_end(td) and lo_end(ta) < lo_end(tb)):
                continue
            Td, Tr = minimal_tree([tc, td]), minimal_tree([ta, tb])
            gd, gr = _gaps(Td, tc, td), _gaps(Tr, ta, tb)
            if any((len(x) > 0) != (len(y) > 0) for x, y in zip(gd, gr)):
                continue
            N = 1 + sum(max(len(x), len(y)) for x, y in zip(gd, gr))
            if best is None or N < best[0]:
                best = (N, tc, td, ta, tb)
    if best is None:
        return None
    N, tc, td, ta, tb = best
    return transport_pair(tc, td, ta, tb), N
