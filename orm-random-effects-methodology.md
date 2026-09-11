# orm Random-Effects Extension: Statistical Methodology

This document is knowledge memory for future development of `rms::orm`'s
random-intercept extension, written for a future session of Claude or
Claude Code picking up this work. It covers the *statistical* side —
what the model is, why it's built this way, and how to interpret and
test it. See `orm-random-effects-implementation.md` for the numerical
algorithms, code architecture, rejected designs, and bug history.

## 1. Overview

`orm.fit`/`orm` fits a single-level random-intercept extension of the
ordinal cumulative-probability (proportional-odds-family) model when
`cluster()` is given, for clustered or repeated-measures data. The base
model, and its `mre`-weighted extension for Markov-type models, are
both described below.

## 2. The base model: single-`sigma` random intercept

For `k` ordinal intercepts and `p` regression coefficients (as in the
ordinary, non-clustered model), a normally-distributed random intercept
`gamma_j ~ N(0, sigma^2)` is added for each distinct cluster (level of
`cluster`). Fitting integrates out the random effects to obtain the
marginal likelihood.

### Why adaptive Gauss-Hermite quadrature (AGQ), not a Laplace approximation

This was empirically driven, not a default choice. Extensive validation
against `ordinal::clmm2` showed that Laplace approximation (equivalent
to `nAGQ=1`) produces substantial, systematic bias in both the fixed
effects and, especially, the variance component — worsening as cluster
informativeness decreases (small clusters, large `sigma`), which is the
regime most common in practice. `nAGQ` itself is not held at a single
fixed value either: it is escalated automatically, refitting
(warm-started) at each successively larger value in a grid
(`c(7,11,15,21,31,45,63)` by default) until the log-likelihood
stabilizes within `nAGQ.tol`. No single fixed `nAGQ` was found adequate
across the full range of `sigma` a user might encounter.

Practical note: the escalation check itself is a *cheap, optional*
verification (it does not drive the fit, only confirms it) — see the
implementation doc for a documented numerical fragility this check can
hit under asymmetric links (`cloglog`), which is benign but worth
knowing about.

## 3. The `mre` extension: weighted random effects (`sigma1`/`sigma2`)

### Motivation

Sometimes the random intercept's influence should not be constant
across every observation in a cluster. The general mechanism: a fixed,
user-supplied per-observation multiplier `mre` scales the random
effect's contribution:

```
eta_ij = ... + [sigma1*(1 - mre_ij) + sigma2*mre_ij] * v_i
```

where `v_i ~ N(0,1)` is a standardized, per-cluster random effect, and
`sigma1`/`sigma2` are two estimated scale parameters. `mre` is entirely
user-computed and fixed — like a column of `x`, never estimated.
`mre=NULL` (the default) recovers the plain single-`sigma` model
exactly.

- `sigma1` is constrained positive via a log link, to pin down `v_i`'s
  otherwise-arbitrary sign (flipping the sign of `sigma2` and every
  `v_i` simultaneously would otherwise leave the likelihood unchanged).
- `sigma2` is left unconstrained — it can legitimately come out
  negative.

### The specific, recommended use case: Markov-1 models

For a Markov-1 model on repeated ordinal responses (a proportional-odds
model whose linear predictor already includes a term for the previous
response, e.g. a spline in the lagged `y`), a natural and recommended
`mre` is a piecewise-constant, `0`/`1` step function: `0` at each
subject's *first* follow-up time, `1` at every subsequent one.

```r
mre <- ave(week, uid, FUN = function(t) as.numeric(t > min(t)))
f <- orm(y ~ x + rcs(prev.y, 4) + cluster(uid), mre=mre, data=d)
```

**Why this matters**: the lagged-response term already captures much
of *why* consecutive observations from the same subject look similar —
that is exactly what an autoregressive term does. A plain,
constant-weight random effect present at *every* visit therefore
double-counts its own influence once the lag term is active: it
re-injects a full-strength correlation contribution at every later
visit on top of what the lag term already explains, which can induce
an unrealistically large or ever-growing correlation across time.
Splitting into `sigma1` (baseline, before the lag term has any prior
value to draw on) and `sigma2` (post-lag) lets the random effect either
add to or subtract from the correlation the lag term alone would
induce, as the data determine — including `sigma2` coming out
negative, which is a legitimate, interpretable correction, not a
symptom of a problem.

`mre` is not restricted to a `0`/`1` step function — any fixed,
user-computed function of the data (continuously varying with
follow-up time, for example) is equally valid.

### Interpreting `sigma1` and `sigma2`

- **`sigma1` large**: substantial subject-level correlation exists at
  baseline, before any lag-term contribution — a real, meaningful
  random effect.
- **`sigma2` small (in magnitude)**, possibly negative: once the lag
  term is doing its job, little *additional* random-effect contribution
  is needed at later visits; a small negative value is the model
  trimming a slight over-induction of correlation that a plain
  constant-weight random effect would otherwise add.

This is exactly the pattern observed in the worked example below.

### Identifiability requirement

`mre` must vary *within* at least some clusters for `sigma1` and
`sigma2` to be separately identifiable. A completely constant `mre`
(all `0`, all `1`, or any other single value across the whole dataset)
leaves them unidentifiable and is rejected with an error. A cluster
contributing only a single observation is harmless — it simply
contributes nothing toward separating `sigma1` from `sigma2`, as long
as *some* clusters elsewhere in the data have within-cluster `mre`
variation.

## 4. Testing: likelihood ratio tests under clustering

### The null model must share the same random-effects structure

A model-level LR test for the covariates (`beta`) requires comparing
the full clustered model against a *clustered* null model — intercepts
plus the *same* random-effects structure (including `mre`, if used),
just without `x`. Comparing against a **non-clustered** null instead
silently tests something else entirely (see the implementation doc's
bug log for the concrete history of this bug) — it does not test
whether the covariates matter, given the random effect is already in
the model.

### The boundary problem in variance-component testing

Under `H0: sigma=0` (or `H0: sigma1=sigma2=0`), the parameter(s) sit at
a *boundary* of their own parameter space, not an interior point. This
is the classic boundary problem in variance-component estimation
(Self & Liang, 1987): the usual asymptotic chi-squared reference
distribution for the LR statistic does not apply as-is, and the
correct reference distribution is a **mixture** of chi-squared
distributions with fewer degrees of freedom than the naive count would
suggest — always making the naive test *conservative*, never
anti-conservative.

Three cases arise in this codebase, each with a different mixture:

| Null hypothesis | Boundary status | Correct reference distribution |
|---|---|---|
| `sigma=0` (plain single-scale model) | 1 parameter, boundary-constrained | 50:50 mixture of chi-sq(0) and chi-sq(1) |
| `sigma1=0` and `sigma2=0`, both boundary-constrained | (Not the actual case here — `sigma2` is unconstrained — included for completeness) | 3-way mixture of chi-sq(0), chi-sq(1), chi-sq(2) |
| `sigma1=0` and `sigma2=0`, only `sigma1` boundary-constrained (the actual `mre` case, since `sigma2` is a regular unconstrained parameter) | 1 boundary + 1 regular parameter | 50:50 mixture of chi-sq(1) and chi-sq(2) |

For the third row (the one that actually applies to the `sigma1`/`sigma2`
model tested against a no-`cluster` null), the critical value at
alpha=0.01 is **8.27**, versus **9.21** for a naive chi-sq(2) — the
correction is real but modest here, since only one of the two
parameters is actually at a boundary.

### Why an MLE at (or near) the boundary and a Bayesian posterior median away from it are not a contradiction

If the marginal likelihood for `sigma` (or `sigma1`) is very flat near
the boundary, a proper prior with support strictly away from zero
(half-Cauchy, half-normal, etc.) can pull the posterior median well
away from zero even when the likelihood's own maximum sits right at
the boundary. This is expected, not a sign either method is wrong — it
reflects how little the *likelihood alone* discriminates among small
values of the variance component in that regime.

## 5. Worked example: the `twstrs` Markov-1 case

This walks through an actual debugging/validation episode, preserved
here because it's a clean illustration of essentially every point
above.

**Setup**: a Markov-1 model (`ptwstrs`, the lagged response, in the
linear predictor) with `cluster(uid)` and no `mre`.

**Symptom**: the single-`sigma` model's `sigma` collapsed cleanly and
monotonically toward the lower boundary across the clustered outer
iterations (`0.0119 -> 8.0e-5 -> 9.8e-8`) — a real MLE-at-the-boundary
result, not an algorithmic artifact (confirmed: the profile search
uses Brent's method on a wide bracket, so this was not a bracket-width
or optimizer-bias issue). Meanwhile, a *Bayesian* fit of essentially
the same model had a posterior median `sigma` of `0.11`.

**Diagnosis**: exactly the mechanism in Section 3 above — `ptwstrs`
already explains most of the within-subject correlation a single,
constantly-weighted random effect would otherwise be needed for, so
the MLE for that one shared `sigma` collapses toward zero. The
Bayesian posterior-vs-MLE discrepancy is exactly the boundary
phenomenon in Section 4.

**Fix**: refit with `mre` distinguishing baseline from post-lag visits:

```r
mre <- ave(week, uid, FUN = function(t) as.numeric(t > min(t)))
f <- orm(twstrs ~ treat*rcs(week,3) + rcs(ptwstrs,4) + rcs(age,4)*sex +
           cluster(uid), mre=mre, data=both, trace=1)
```

**Result**: `sigma1` came out large, `sigma2` around `-0.2` — exactly
the interpretation in Section 3 (substantial baseline correlation,
small negative post-lag correction) — and the clustered full model
beat the *no-cluster-at-all* model by a deviance of `10`. Since
`sigma1=sigma2=0` (no cluster term at all) is the one-boundary/
one-regular-parameter case from Section 4's table, the correct
critical value at alpha=0.01 is `8.27` — `10` clears it comfortably
(`p ≈ 0.004` under the correct mixture, `p ≈ 0.007` under a naive
chi-sq(2)), so the improvement is real, not an artifact of adding two
free parameters or of using the wrong reference distribution.

## 6. Flagged for future examination: interaction with y-dependent effects (yde)

**This has not yet been examined and needs dedicated attention before
combining with `cluster()`/`mre`.** `orm`/`lrm` support *y-dependent
effects* (also called partial or constrained partial proportional-odds
models): certain covariates are allowed to have effects that differ by
which side of a response-category cutpoint an observation falls on,
i.e. they violate the proportional-odds assumption in a controlled way.
The underlying Fortran interface (`ormll`/`ormeta`) already carries the
machinery for this — separate `lp1`/`lp2` linear predictors and
separate `ia`/`ia2` index arrays keyed to which side of a cutpoint an
observation's response falls on — but **every current call site in the
random-effects code (`orm.rfit.r`) sets `lp1=lp2=lp` unconditionally**.
The entire `cluster()`/`mre`/AGQ machinery has only ever been built,
tested, and validated under a full proportional-odds assumption.

Before allowing y-dependent effects together with `cluster()`, a future
session needs to examine (see the implementation doc, Section 8, for
the specific technical questions this raises in `clusterModeFind`,
`agqStep`, and `sparseMissingInfo`).

## References

- Self, S.G. and Liang, K-Y. (1987). Asymptotic properties of maximum
  likelihood estimators and likelihood ratio tests under nonstandard
  conditions. *JASA* 82(398), 605-610. — the boundary-testing mixture
  chi-squared results used in Section 4.
- `ordinal::clmm2` — used as the empirical reference for validating
  AGQ vs. Laplace bias (Section 2).
