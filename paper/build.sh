#!/bin/sh
# Build paper/main.pdf (pdflatex only; the bibliography is inline).
set -e
cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
