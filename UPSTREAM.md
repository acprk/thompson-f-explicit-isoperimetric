# Upstream provenance

This repository contains no code or text copied from other repositories. It analyses one
preprint, and its Lean formalization, in [OpenAI's mathematics repository](https://github.com/openai/math),
at commit
[`fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb`](https://github.com/openai/math/tree/fd4aeeb2ee4fc729c18d98444fed42fd0529eeeb):

- *Thompson's group F is nonamenable* (OpenAI, September 23, 2026), directory
  `preprints/Thompsons-group-F-is-nonamenable-September-23-2026`.
- Lean formalization, directory `lean/OAI/GroupTheory/Thompson/`. The note refers to these files:
  - `lean/OAI/GroupTheory/Thompson/AnalyticNormalization.lean`
    (`normalizeNeg_lipschitz`, `normalizedMap_lipschitz`: the `2/ε` normalization bound);
  - `lean/OAI/GroupTheory/Thompson/AnalyticConstruction.lean`
    (the separation constant obtained by `Classical.choose`, and `displacementLipschitzConstant`);
  - `lean/OAI/GroupTheory/Thompson/TubeRegularCurve.lean` (`RegularCurveData`, tube radius
    `ρ = min(1/1024, σ/8)`);
  - `lean/OAI/GroupTheory/Thompson/AnalyticCorrection.lean` (`CorrectionData.lipschitzConstant`);
  - `lean/OAI/GroupTheory/Thompson/CorrectionLocal.lean`.

All of the mathematics of nonamenability is in that preprint. This repository only tracks its
constants and makes them explicit for the generating set `{x0, x1}`. We have not modified the
Lean code. Section and statement numbers quoted in the note (Sections 2.1–2.3, Appendix A)
refer to the preprint at the commit above.
