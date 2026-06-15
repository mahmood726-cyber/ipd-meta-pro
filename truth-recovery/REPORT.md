# Truth-Recovery Validation -- ipd-meta-pro

**Verdict: SHIP (genuine IPD meta-analysis engine; both one-stage and two-stage
recover truth, with one honest, well-understood coverage caveat for one-stage).**

## What this repo is

A large modular single-file IPD meta-analysis app. The statistical core lives in
`dev/modules/`. Two genuine estimators of an IPD treatment effect are present and
were validated against a known simulation truth:

- **One-stage**: `OneStageGLMM.fitBinary` (`02_25_one-stage-glmm.js`) -- a
  penalised quasi-likelihood (Breslow-Clayton) logistic GLMM with **random study
  intercepts** and treatment as a fixed effect.
- **Two-stage**: per-study log-OR (2x2 -> Woolf, 0.5 correction only on zero
  cells) pooled with `MetaAnalysis.randomEffectsREML` and
  `MetaAnalysis.applyHKSJ` (`02_08_meta-analysis.js`), the latter using
  t_{k-1} (Hartung-Knapp-Sidik-Jonkman).

## Method

`truth-recovery/dgp-ipd.mjs` is a STANDALONE seeded DGP (independent RNG) using
the canonical random-effects logistic generative model:
study baseline `alpha_i = alpha0 + N(0,sigmaAlpha^2)`, study treatment effect
`theta_i = theta + N(0,tau^2)`, `logit P(Y=1) = alpha_i + theta_i*t`. The KNOWN
TRUTH is the overall treatment log-OR `theta` and between-study SD `tau`;
individual binary outcomes are generated.

`truth-recovery/harness.mjs` runs the repo's OWN one-stage and two-stage fits
(imported verbatim from `engine.mjs`, which bundles modules 02_06, 02_07, 02_08,
02_25 unedited) over 200 Monte-Carlo replicates per scenario and measures bias of
the recovered log-OR plus empirical coverage of the true log-OR.

## Results (200 replicates/scenario; cov = coverage of true log-OR)

| scenario              | method         | mean est | bias   | coverage | mean CI width |
|-----------------------|----------------|---------:|-------:|---------:|--------------:|
| k=10, tau=0.1 (OR0.7) | one-stage GLMM | -0.357   | -0.001 | 91.5%    | 0.24 |
|                       | two-stage z    | -0.358   | -0.001 | 96.0%    | 0.28 |
|                       | two-stage HKSJ | -0.358   | -0.001 | 99.0%    | 0.33 |
| k=10, tau=0.3 (null)  | one-stage GLMM | +0.007   | +0.007 | 73.5%    | 0.23 |
|                       | two-stage z    | -0.004   | -0.004 | 95.0%    | 0.42 |
|                       | two-stage HKSJ | -0.004   | -0.004 | 96.5%    | 0.49 |
| k=5,  tau=0.3 (OR0.7) | one-stage GLMM | -0.332   | +0.025 | 80.5%    | 0.42 |
|                       | two-stage z    | -0.343   | +0.013 | 90.5%    | 0.66 |
|                       | two-stage HKSJ | -0.343   | +0.013 | 99.5%    | 1.18 |

## Findings (classic IPD-MA questions, as measured)

**(a) One-stage vs two-stage agreement.** Point estimates agree closely
everywhere (mean-estimate gap < 0.015 log-OR). Both are essentially unbiased for
the true log-OR; small-k shrinks bias to ~0.013-0.025 (expected finite-sample
noise). The two methods *recover the same truth* in their point estimate.

**(b) Clustering / random-effects structure drives the SE, and the one-stage
model under-covers under heterogeneity.** The repo's one-stage GLMM (`fitPQL`)
fits **random study intercepts only** -- treatment is a fixed effect with no
random slope. When the truth has genuine between-study treatment-effect
heterogeneity (tau>0), the GLMM's SE reflects only the intercept variance and is
therefore too small, so coverage drops below nominal: 91.5% at tau=0.1, and
**73.5% (k=10) / 80.5% (k=5) at tau=0.3**. This is the well-documented PQL /
random-intercept-only limitation, reproduced here from the repo's own code -- not
a coding bug, but a model-specification caveat users must know. The two-stage
random-effects model, which propagates tau^2 into the pooled SE, is correctly
calibrated at k=10 (95-96%).

**(c) Small-k coverage: HKSJ vs naive z.** At k=5 with heterogeneity, two-stage
naive z under-covers (90.5%); the repo's `applyHKSJ` (t_{k-1}) repairs this,
raising coverage to 99.5% -- conservative/over-covering, as HKSJ is known to be
when Q is small, but on the safe side. This is the textbook small-k result and
the repo implements it correctly.

**Which approach recovers truth best?** For the *point estimate*, both are
equivalent. For *valid inference* (CI that covers the truth), the **two-stage
random-effects model -- with HKSJ at small k -- is the safe default** in these
scenarios. The one-stage GLMM is competitive only when between-study
treatment-effect heterogeneity is negligible; with tau>0 it is anticonservative
because it lacks a random treatment slope.

## Recommendation

**SHIP.** Both estimators are genuine and recover the known true log-OR. Validation
is additive (`truth-recovery/` only); no source module edited. Recommend a
user-facing note (and a future enhancement) that the one-stage GLMM uses random
intercepts only and can under-cover when treatment-effect heterogeneity is
present; prefer the two-stage RE(+HKSJ) result, or add a random treatment slope to
the one-stage model, in that regime.

## Reproduce

    node truth-recovery/harness.mjs                       # results table
    node --test truth-recovery/test-truth-recovery.mjs    # 7 assertions
