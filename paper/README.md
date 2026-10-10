# Paper

`main.tex` is the note (11 pages); `main.pdf` is the compiled version. Build with `./build.sh`
(pdflatex, run twice; the bibliography is inline and needs no biber). It uses `pgfplots` for
Figure 1. The author block is a placeholder.

Scripts named in the note live in [`../computations/`](../computations/): `cert.py`,
`s1_sigma_paper_curve.py` and `s4_wordlength_R.py` at the top level, and `v1_sigma.py`,
`v2_words.py`, `v3_avg.py`, `v4_constants.py` and `r1_F2sharp.py` in `computations/independent/`.
The upstream preprint and Lean files cited as [OAI] are listed in [`../UPSTREAM.md`](../UPSTREAM.md).
