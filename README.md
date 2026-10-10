# An explicit isoperimetric constant for Thompson's group F

OpenAI's preprint *Thompson's group F is nonamenable* (September 23, 2026) proves
nonamenability by a finite averaging argument, and formalizes it in Lean. Neither the preprint
nor the formalization gives an explicit boundary constant, or a constant for a prescribed
generating set.

This repository contains a short note and the computations behind it. For the standard
generators `x0, x1` and the edge-isoperimetric (Cheeger) constant of the 4-regular Cayley graph

```
h_E(F) = inf_{A ⊂ F finite, nonempty} ( |x0 A △ A| + |x1 A △ A| ) / |A| ,
```

it proves:

- **unconditionally**, `h_E(F) ≥ 2.109 · 10⁻¹⁹`. As far as we know this is the first explicit
  positive lower bound for an isoperimetric constant of `F`;
- with a sharper word-length lemma (F2′), `h_E(F) ≥ 2.466 · 10⁻¹⁹`;
- **conditionally**: if the closed unit ball of an infinite-dimensional Hilbert space admits a
  `k`-Lipschitz retraction onto its sphere, then `h_E(F) ≥ g(k)` for an explicit function `g`.
  At `k = 29` this gives `h_E(F) ≥ 8.596 · 10⁻⁹`. The bound `k₀(ℓ²) ≤ 28.99` is quoted in the
  literature from Baronti–Casini–Franchetti (2003) and Casini–Piasecki (2017). **We have read
  neither paper; we rely on second-hand quotations.**

The best known upper bound is `h_E(F) < 0.4996` (Guba's density bound `> 3.5004`).

The results were found and checked by **Claude (Anthropic) agents** working under the direction
of the repository owner, [AUTHOR NAME]. **This is a standard effectivization of the upstream
argument, not new mathematics.** It is conditional on the correctness of the upstream preprint.
No human peer review has taken place, and nothing here is formalized in Lean. See
[Status and caveats](#status-and-caveats).

**Note:** [*An explicit isoperimetric constant for Thompson's group F*](paper/main.pdf)
([LaTeX source](paper/main.tex), [build instructions](paper/README.md)).

## What is proved

| Statement | Bound | Kind of proof |
| --- | --- | --- |
| `F` is nonamenable, i.e. `h_E(F) > 0` | no explicit constant | upstream (preprint + Lean) |
| Separation constant of the upstream curve: `σ(1/32) = 1/256` (**Prop. 4.1**) | exact | note; computer-assisted (Arb branch-and-bound), not formalized |
| Hilbert-space normalization is `1/ε`-Lipschitz, not `2/ε` (**Lemma 2.2**) | factor 2 in `L` | note, by hand (two-line identity) |
| `h_E(F) ≥ c⋆` (**Theorem 1.1**), unconditional | `c⋆ = 2.1091… · 10⁻¹⁹` | note, by hand; constants certified in exact rationals / Arb |
| Same with the sharper word count (F2′, Lemma 5.3) (**Theorem 6.1**) | `2.4666… · 10⁻¹⁹` | note, by hand; (F2′) checked by computer |
| If a `k`-Lipschitz retraction `B → S` exists, `h_E(F) ≥ g(k)` (**Theorem 1.2**) | `g(29) = 8.5964 · 10⁻⁹` (`g′(29) = 1.0197 · 10⁻⁸`) | note, by hand; conditional |
| Spectral form (**Corollary 1.3**): `λ₀ ≥ h_E²/8 ≥ 5.56 · 10⁻³⁹` | | Mohar's inequality + Theorem 1.1 |
| Ceiling of this method: `h_E ≥ 2.51 · 10⁻⁴` at best | | note (Goebel's displacement bound) |

For comparison, keeping the `2/ε` normalization constant of the Lean code gives
`5.10 · 10⁻²⁰`; the `1/ε` constant improves this by a factor `4.13`.

### Where the constant was lost

Ledger for Theorem 1.1, in binary orders (`log₂`):

| Factor | `log₂` | Origin |
| --- | --- | --- |
| `δ² = 1/4` | −2.00 | displacement of `f = −G/‖G‖` |
| `ε² = (ρ/4)²`, `ρ = 1/2048` | −26.00 | `ρ = σ/8`, `σ = 1/256 = (1/8)·(1/32)` (curve speed × separation scale) |
| `Lip(G)⁻²`, `Lip(G) = 3615.86` | −23.64 | dominated by the cutoff term `7/(4ρ) = 3584` |
| **`κ = (δ/L)²`** | **−51.64** | analytic map |
| `(1 − 4/(Dκ)) / (4(1 − 1/D))`, `D = 2⁵⁹` | −2.04 | variance constant 4 |
| `1/M̄`, `M̄ = 329.5` | −8.36 | average length of the transport words, about `5.5m` |
| **`h_E ≥ c⋆`** | **−62.04** | (−61.81 with (F2′)) |

About 52 of the 62 binary orders come from the analytic displacement map, through
`L ≈ 7/ρ²`. The group-theoretic part, which is the only part that depends on `F`, costs
about 10, of which 8.4 are the average word length of the transports in `x0, x1`.

The note makes three routine changes to the upstream argument, and nothing else:

1. it evaluates the one non-explicit constant of the Lean construction (a curve-separation
   constant obtained by `Classical.choose`) as `σ(1/32) = 1/256`;
2. it uses the Hilbert-space Lipschitz constant `1/ε` of radial normalization instead of the
   Banach-space constant `2/ε`;
3. it replaces the finite set `S` of affine transports by `{x0, x1}`. It uses explicit
   Cannon–Floyd–Parry normal-form words, and averages the transport errors instead of taking
   their maximum.

## Repository layout

- `paper/`: LaTeX source and PDF of the note (11 pages), plus `build.sh`.
- `computations/cert.py`: the main certificate (about 15 s, one core). Every number that enters
  a stated constant is computed in exact rational arithmetic or in Arb ball arithmetic
  (`python-flint`, outward rounding). It covers:
  - [A] `σ(1/32) = 1/256` (analytic cases and a 6462-box branch-and-bound);
  - [B] the constant chain `ρ`, `Lip(G)`, `L`, `κ`;
  - [C] the normal form, (F1), (F2), (F2′) and the closed-form averages: exhaustive for
    `m ≤ 4`, random pairs up to `m = 64`, including the transport property of the words;
  - [D] all bounds of Theorems 1.1, 1.2 and 6.1 and of Table 1;
  - [E] the ceilings, the ledger and the spectral corollary.

  Its output is in `computations/expected_output/cert_output.txt`.
- `computations/`: further scripts written separately from `cert.py`:
  - `thompson.py`: exact toolkit for `F`;
  - `s1_sigma_paper_curve.py`: a second Arb proof of `σ(1/32) = 1/256`;
  - `s2_constants.py`: the constants with the `2/ε` normalization;
  - `s4_wordlength_R.py`: (F1), (F2) and the closed forms, exhaustive for small `m`;
  - `s3_explore.py`: a **non-rigorous Monte Carlo** exploration of alternative transports,
    which the note does not use.
- `computations/independent/`: re-implementations with their own data structures:
  - `v1_sigma.py`: `σ(1/32) ≥ 1/256` in `mpmath.iv`, with a different region split;
  - `v2_words.py`, `v2_test.py`, `v2b_test.py`: (F1), (F2) and the transport property, on
    random pairs up to `m = 64`;
  - `v3_avg.py`: the closed-form averages against exhaustive enumeration;
  - `v4_constants.py`: the (F2)-based constants;
  - `r1_F2sharp.py`: (F1), (F2′), the averages, and whole-PL-map transport checks.
- `computations/run_all.sh`: runs `cert.py` (and compares its output with the saved one), then
  the fast checks, in about 20 s. `FULL=1 ./run_all.sh` adds the slow word-length checks
  (about 1 h). Requires Python 3 with `python-flint` (0.9) and `mpmath`.

## Status and caveats

- **Standard effectivization; no new mechanism.** No step of the upstream argument fails or needs
  a new idea. The only contribution is that the constant is explicit.
- **Computer-assisted step.** `σ(1/32) = 1/256` (Proposition 4.1) is proved by a certified
  interval computation, checked by three separate implementations. It is not formalized.
  Substituting it for the `Classical.choose` value in the Lean development is justified in the
  note, but we have not changed or re-run the Lean code.
- **The conditional bound rests on a source we have not read.** We have not seen the full text
  of Baronti–Casini–Franchetti (*J. Approx. Theory* 120, 2003) or of Casini–Piasecki
  (*J. Nonlinear Convex Anal.* 18, 2017). We rely on the quotations by
  Chaoha–Goebel–Termwuttipong and by Mahesh, so we have not checked whether `28.99` is a proved
  theorem or a numerical estimate. We state the `8.596 · 10⁻⁹` bound only conditionally.
- **References.** Volume, page and year data were checked against publisher pages for Guba
  (2004, 2025), Belk–Brown, Mohar, Chaoha–Goebel–Termwuttipong and Baronti–Casini–Franchetti;
  the statements of Mohar (1988) and of the retraction-constant papers were not read in full.
- **Formalization scope.** Only the analytic displacement map is taken from the upstream Lean
  development; the dyadic layout, the averaged form of the counting proposition and the
  word-length lemmas are proved in the note and checked by computer, not in Lean.
- **Conditional on upstream.** We use the upstream preprint and its Lean formalization as given.
  If the upstream proof has a gap, so do our theorems.
- **No human peer review. Not formalized in Lean.**

## AI attribution

The results in this repository were found, verified and written up by Claude (Anthropic) agents
working under the direction of the repository owner, [AUTHOR NAME]. This covers the analysis,
the note and the computations, including the separately written check scripts. The human owner
directed the work and is responsible for its release.

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite the OpenAI preprint listed in
[`UPSTREAM.md`](UPSTREAM.md).
