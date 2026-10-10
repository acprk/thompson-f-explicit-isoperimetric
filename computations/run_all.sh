#!/bin/sh
# Reproduce the computations.
# Requires Python 3 with python-flint (0.9, Arb ball arithmetic) and mpmath:
#     pip install python-flint mpmath
# Default: the main certificate and the fast cross-checks (about 20 s on one core).
# FULL=1 ./run_all.sh additionally runs the slow exhaustive/random word-length checks
# (roughly 1 h in total; r1_F2sharp.py uses 4 processes).
set -e
cd "$(dirname "$0")"

echo "== cert.py: main certificate, all stated constants (~15 s)"
python3 cert.py > /tmp/cert_output.$$.txt
if diff -q /tmp/cert_output.$$.txt expected_output/cert_output.txt >/dev/null; then
  echo "   output identical to expected_output/cert_output.txt"
else
  echo "   WARNING: output differs from expected_output/cert_output.txt"; exit 1
fi
rm -f /tmp/cert_output.$$.txt

echo "== thompson.py: self-test of the tree-pair / normal-form toolkit (<1 s)"; python3 thompson.py
echo "== s1_sigma_paper_curve.py: sigma(1/32) = 1/256, Arb branch-and-bound (<1 s)"; python3 s1_sigma_paper_curve.py | tail -3
echo "== s2_constants.py: constants with the 2/eps normalization (<1 s)"; python3 s2_constants.py | tail -5
echo "== independent/v1_sigma.py: sigma(1/32) >= 1/256 in mpmath.iv (~2 s)"; (cd independent && python3 v1_sigma.py)
echo "== independent/v3_avg.py: closed-form averages vs exhaustive enumeration, m <= 4 (~2 s)"; (cd independent && python3 v3_avg.py)
echo "== independent/v4_constants.py: (F2)-based constants (<1 s)"; (cd independent && python3 v4_constants.py)

if [ "${FULL:-0}" = 1 ]; then
  echo "== s4_wordlength_R.py: (F1), (F2), transports, exhaustive m <= 4 (~3 min)"; python3 s4_wordlength_R.py | tail -4
  echo "== independent/v2_test.py: (F1), (F2), transports, random pairs (~7 min)"; (cd independent && python3 v2_test.py)
  echo "== independent/v2b_test.py: (F1), (F2), transports, random pairs up to m=64 (~9 min)"; (cd independent && python3 v2b_test.py)
  echo "== independent/r1_F2sharp.py: (F1), (F2'), averages, full-map transports (~35 min wall, 4 processes)"; (cd independent && python3 r1_F2sharp.py)
  echo "== s3_explore.py: NON-RIGOROUS Monte Carlo exploration, not used in the note (~2.5 min)"; python3 s3_explore.py
fi
echo "done"
