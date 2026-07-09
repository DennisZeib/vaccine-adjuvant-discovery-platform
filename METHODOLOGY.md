# Methodology & Limitations

## What this is

A computational simulation and scoring tool for exploring microbial
co-culture dynamics and their relevance to candidate vaccine adjuvants
(PAMP targets: LPS, mannan, flagellin). It combines:

- An ODE-based kinetic model (Monod growth + Luedeking-Piret production)
  simulating microbial co-culture over time.
- A Hill-equation dose-response scorer estimating adjuvant potential
  from simulated metabolite concentration.
- Monte Carlo sensitivity analysis and Pareto-front optimization across
  simulation parameters (glucose, duration, inoculation ratio).

## What this is NOT

This is **not** a validated pharmacological or immunological prediction
tool. It cannot tell you whether a real adjuvant candidate is safe or
effective in vivo. Specifically, it does not model:

- Immune memory or adaptive response kinetics beyond a single-dose
  static scoring heuristic.
- Route of administration, formulation chemistry, or
  adjuvant-antigen co-formulation effects.
- Species-specific immune variation, age, comorbidities, or genetic
  background.
- Correlates of protection — i.e. whether a high "score" here
  correlates with real-world protective immunity.
- Real pharmacokinetics/pharmacodynamics (absorption, distribution,
  metabolism, excretion).

Any output from this tool should be treated as a **hypothesis-generation
aid at best** — a way to explore parameter sensitivity and relative
tradeoffs, not a substitute for in vitro assays, animal studies, or
clinical evaluation.

## Units: the most important caveat

The simulator's `production` coefficients (alpha/beta, Luedeking-Piret
model) and resulting metabolite concentrations are expressed in
**simulator-internal "toy units"** — they were tuned so the ODE system
produces numerically stable, order-of-magnitude-plausible dynamics, not
because they map onto any real measured unit (mM, μg/mL, etc.).

This matters because real literature dose-response data (e.g. EC50
values reported in picomolar or microgram doses) is **not directly
compatible** with this scale. During literature review (July 2026), we
found:

- Flagellin: real TLR5-binding EC50 = 2.4 ± 1.4 pM (Smith et al. 2013,
  DOI: 10.1002/bit.24903); human vaccine dose-response threshold
  around 6-10 mcg (DOI: 10.1016/j.vaccine.2017.09.070).
- Mannan-MUC1: mouse dose-response showing 1-7 mcg favoring cellular
  immunity, >7 mcg favoring humoral immunity (DOI: 10.1007/s002620050449).
- LPS: no real TLR4-binding potency data was located during search;
  most literature addresses LPS as a stimulus for testing *other*
  compounds, not LPS's own dose-response curve.

These real values are recorded as `citation` fields in `_KNOWLEDGE_BASE`
for reference and future recalibration work, but the `ec50` field used
by the scorer remains simulator-internal until a full unit-system
recalibration is done (see "Future work" below). **Do not interpret
`ec50` values in this codebase as real-world potency figures.**

## Known bug history (for transparency)

An earlier version of this tool had a bug where Monte Carlo sensitivity
analysis returned `std ≈ 0.0` across all samples — i.e. varying inputs
(glucose, duration, ratio) had no effect on the output score. Root
cause: `ec50` values (5.0, 3.0, 8.0 in original toy units) were far
below the concentrations the simulator actually produced (roughly
80-100, 22-44, and 4-16 respectively for LPS, flagellin, and mannan),
causing the Hill dose-response curve to saturate near 1.0 regardless of
input variation. Fixed by recalibrating `ec50` to the simulator's
actual observed concentration ranges (commit `d4a6bac` and prior).

## Future work

- Full unit-system recalibration: define what one "toy unit" of
  simulated concentration is meant to represent in real terms (e.g.
  mg/mL equivalent), and rescale either the kinetic model or the
  dose-response layer so real literature EC50 values can be used
  directly, rather than as documentation-only citations.
- Literature-source EC50/kinetics for LPS specifically — currently
  the largest gap.
- Validate `mu_max`/`Ks`/`Y_xs` growth parameters against published
  growth-curve data for at least the most commonly used organisms
  (E. coli K-12 is the best-studied candidate for this).
- Add explicit uncertainty ranges (not just point estimates) once
  literature EC50s are properly integrated, reflecting real
  study-to-study variability.

## Intended audience & use

Built as an educational/exploratory tool and portfolio project
demonstrating computational modeling, Monte Carlo methods, and
multi-objective optimization. Feedback from anyone with immunology,
vaccinology, or computational biology expertise is welcome — please
open an issue or reach out directly.
