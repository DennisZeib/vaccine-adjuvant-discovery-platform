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

## Growth kinetics correction: S. cerevisiae Y_xs (verified by simulation)

Unlike E. coli K-12 (where the search only confirmed plausibility), a
targeted search for S. cerevisiae growth kinetics found a real,
specific discrepancy: the simulator's `Y_xs = 0.40` was roughly 2-4x
higher than published wild-type biomass yield data (0.11 g/g glucose
under Crabtree-positive conditions, up to 0.21 g/g in engineered
Crabtree-relieved strains — DOI: 10.1016/j.synbio.2022.06.004).

**Fix applied:** `Y_xs` recalibrated from 0.40 to 0.20 (documented
middle estimate, not a precise fit — see inline citation in
`create_default_microbes()`).

**Verified by re-running the simulation** (not just a code change
left unverified): Monte Carlo comparison of Saccharomyces cerevisiae
+ Bacillus licheniformis (glucose=100, duration=48h, n=60) after the
fix showed mean score 61.4, std 1.0, with toxicity values in the
79-126 range — consistent with the expected direction of the change
(lower biomass yield per glucose unit shifts more carbon toward
overflow/byproduct metabolites). Parameter sensitivity analysis on
this run showed vaccination ratio as the dominant driver (Pearson
r = -0.80), more than 3x the effect of culture duration (r = +0.24),
a specific and non-obvious result that would not have been trustworthy
before both the ec50 saturation bug and this yield-coefficient
placeholder were corrected.

## Growth kinetics validation (E. coli K-12)

A targeted EuropePMC search (July 2026) checked the simulator's default
E. coli K-12 growth parameters against published chemostat/continuous-
culture literature:

**Current simulator values** (`create_default_microbes()`):
`mu_max = 0.85 /h`, `Ks = 5.0 mM`, `Y_xs = 0.50`, `k_death = 0.015 /h`.

**mu_max — plausible, weakly supported.** A 2005 glucose-limited
chemostat study (E. coli K-12 strain TG1, 28°C) observed washout
occurring near a dilution rate of 0.4/h (DOI: 10.1099/mic.0.27481-0).
Since a chemostat's dilution rate cannot exceed `mu_max` without
washout, this implies a real `mu_max` somewhat above 0.4/h under those
(cooler, minimal-media) conditions — consistent with commonly cited
37°C glucose `mu_max` values in the 0.7-1.0/h range reported elsewhere
in the microbiology literature. The simulator's 0.85/h sits within this
plausible range.

**Y_xs — not contradicted, not confirmed.** No single clean yield-
coefficient value was found in this search batch; multiple fed-batch
and chemostat studies discuss yield coefficients in a broadly similar
range without providing a directly comparable reference number.

**Ks — unresolved, flagged uncertain.** No glucose half-saturation
constant was found in this literature batch. Other published estimates
for E. coli's glucose Ks (outside this specific search) are commonly
cited in the 0.01-0.2 mM range, well below the simulator's `Ks = 5.0`.
Given the simulator's concentrations are toy units rather than real mM
(see Units section above), this discrepancy may not have the same
practical impact as the LPS ec50 mismatch did, but it has not been
resolved and should not be assumed correct.

**k_death** was not addressed by this literature batch.

**Conclusion:** growth kinetics were built with plausible order-of-
magnitude values for the best-supported parameter (`mu_max`), unlike
the arbitrary placeholder `ec50` values that caused the earlier
Monte Carlo bug. `Ks` remains the most uncertain parameter and is a
candidate for future targeted literature search.

## Future work

- Full unit-system recalibration: define what one "toy unit" of
  simulated concentration is meant to represent in real terms (e.g.
  mg/mL equivalent), and rescale either the kinetic model or the
  dose-response layer so real literature EC50 values can be used
  directly, rather than as documentation-only citations.
- Literature-source EC50/kinetics for LPS specifically — currently
  the largest gap (see "LPS EC50 search" below for what was tried).
- Resolve E. coli K-12's `Ks` value against dedicated glucose
  half-saturation literature (see "Growth kinetics validation" above).
- Validate growth kinetics for the other 15 microbes in
  `create_default_microbes()` — only E. coli K-12 has been checked
  so far.
- Add explicit uncertainty ranges (not just point estimates) once
  literature EC50s are properly integrated, reflecting real
  study-to-study variability.

## LPS EC50 search: a documented gap, not a skipped step

Two rounds of targeted EuropePMC literature search (July 2026) failed
to locate a clean, absolute EC50 value (in molar or mass-concentration
units) for native LPS's own TLR4/MD-2 binding potency — despite finding
such values readily for flagellin (TLR5) and mannan-MUC1 (see
`_KNOWLEDGE_BASE` citations above).

**What was searched:** general LPS+TLR4+EC50 queries, then narrowed to
the standardized E. coli O111:B4 reference material specifically (the
most widely used reference endotoxin in immunology, with a WHO/NIBSC
reference standard).

**What was found instead:**
- Relative sensitivity data: TLR4 receptor overexpression shifts LPS's
  effective EC50 30-fold left; the LPS-resistant Tlr4(Lps-d) mutant
  (found in C3H/HeJ mice) shifts it 2600-fold right (Poltorak et al.
  1999, DOI: 10.1006/bcmd.1999.0262). These are fold-shifts under
  genetic manipulation, not a baseline absolute value.
- A 1983 paper (DOI: 10.1016/s0092-1157(83)80013-4) establishing a
  *reference endotoxin* potency standard via collaborative Limulus
  amebocyte lysate (LAL) gelation assay across 14 laboratories — this
  confirms LPS potency is standardized in the field, but as an
  operational/comparative gelation-endpoint unit, not a molar EC50.
- A validated methodology exists (TLR4 HEK reporter-gene EC50 assays)
  for measuring this directly, suggesting the number likely exists in
  papers this search did not surface — this remains open for future
  literature work.

**Working hypothesis for why this is harder than flagellin:** LPS
structure (particularly the lipid A moiety) varies substantially
across Gram-negative bacterial species and even strains, unlike
flagellin's more evolutionarily conserved TLR5-recognition domain.
This may make a single universal "LPS EC50" less standard in the
literature than a strain-specific one. This hypothesis is not yet
verified against dedicated structural biology literature.

**Current state:** `_KNOWLEDGE_BASE["polysaccharide_lps"]["ec50"]`
remains simulator-internal (see Units section above), with no real
citation available as of this writing. This is flagged deliberately
rather than silently left unlabeled.

## Intended audience & use

Built as an educational/exploratory tool and portfolio project
demonstrating computational modeling, Monte Carlo methods, and
multi-objective optimization. Feedback from anyone with immunology,
vaccinology, or computational biology expertise is welcome — please
open an issue or reach out directly.