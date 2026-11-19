# AI Use Statement

This document records how generative AI and AI-assisted tools were used in
this repository. It complements the declaration in the companion manuscript,
*Simulation and maximum likelihood estimation for a stochastic
convection-diffusion equation*, and is organized after the **Diligence**
competency of the 4D AI Fluency framework: creation, transparency, and
deployment diligence.

## Tools

| Tool | Scope of use | Trace in this repository |
| --- | --- | --- |
| Warp terminal, Oz agent | Doxygen and GitHub Actions documentation workflows; snapshot plotting and 3-D rugosity targets | `Co-Authored-By: Oz` trailer on `a1f1265`, `57e53c0`, `0275fc0`, `bee2d99` |
| Gemini Code Assist | Manuscript-quality figure styling in `visualization/scripts/` | `.aiexclude` |
| OpenAI Codex | Software engineering tasks | `.codex/` |
| Grammarly | Grammar, spelling, and readability of the manuscript text only; no scientific content | Manuscript, not this repository |
| Claude Code (Anthropic) | Drafting this statement and the manuscript AI declaration; the skew-symmetry unit test; license and citation metadata | This file; `tests/unit/test_core.f90`; `LICENSE`, `CITATION.cff`, `.zenodo.json` |
| Unattributed LLM sessions (OpenAI Codex) | Progress/ETA reporting and run logging; ParaView noise-effect video | `IMPLEMENTATION_SUMMARY.md`, `SESSION_RESUME_*` |

## Creation diligence: what AI was and was not used for

AI assistance was limited to software engineering: refactoring, progress
reporting, logging, build and CI configuration, documentation, and
post-processing and visualization scripts. The only exception is
Grammarly, which was used for grammar, spelling, and readability of the
manuscript text.

AI was **not** used to formulate the model, derive the spectral Galerkin
system, derive the exact entries of the advection matrix, derive the maximum
likelihood estimators, or prove their consistency and asymptotic normality.
The numerical kernels in `src/modules/spectral_operators_mod.f90`,
`sde_simulation_mod.f90`, and `parameter_ml_estimation_mod.f90` implement
formulas derived by the authors.

Large data, logs, and generated outputs are excluded from AI indexing through
`.aiexclude`.

## Transparency diligence: how AI contributions are disclosed

- Commits produced with an agent should carry a `Co-Authored-By:` trailer
  naming the tool.
- Work done in an AI session that is not committed directly, such as figure
  iterations or video renders, is summarized in a `SESSION_RESUME_*` file or
  in the commit message that lands it.
- AI tools are never listed as authors of the software or the manuscript.

Commits made before this policy existed are only partially attributed. The
table above is the authoritative record for them.

## Deployment diligence: how AI output is verified

The authors reviewed all AI-generated code before merging it, and remain
responsible for it. In addition:

1. `make test` (unit and end-to-end smoke tests) must pass locally and in the
   GitHub Actions `CI` workflow.
2. Structural properties that the theory depends on are tested directly. For
   example, `test_interaction_matrix_skew_symmetry` checks that the advection
   matrix satisfies `max|A + Aᵀ| / max|A| ≤ 64 ε` for square and rectangular
   domains and for several velocity modes. On the 50 × 50 configuration used
   in the manuscript the measured defect is about `1e-16`.
3. Visualization changes do not modify simulation data. Every figure and video
   can be regenerated from the archived data with the scripts in
   `visualization/scripts/`.
4. The results reported in the manuscript come from the reviewed code archived
   at [DOI 10.5281/zenodo.21252400](https://doi.org/10.5281/zenodo.21252400).
