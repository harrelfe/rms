# orm Random-Effects Extension: Implementation Notes

This document is knowledge memory for future development of `rms::orm`'s
random-intercept extension, written for a future session of Claude or
Claude Code picking up this work. It covers the *numerical/code* side —
architecture, algorithms tried and rejected, and a bug history with
enough detail to avoid re-discovering the same failure modes. See
`orm-random-effects-methodology.md` for the statistical model, its
motivation, and how to interpret and test it.

**Key files**: `orm.fit.r` (dispatch: fixed-effects fit, cluster
dispatch, LR-test/deviance assembly), `orm.rfit.r` (the random-effects
fitting workhorse: `clusterModeFind`, `agqStep`, `ormrfit`,
`sparseInfoMatrix`, `sparseMissingInfo`), `infoMxop.r` (information
matrix assembly/inversion, including the `sigma1`/`sigma2` covariance
extraction). The compiled routines `ormll`/`ormeta`/`ormidx` live in
Fortran (`ormll.f90` and a separate `ormeta.f90`) — **treat any copy of
these `.f90` files or `init.c` in a snapshot as potentially stale**;
this was confirmed to happen mid-project (the developer works across
two machines, and `init.c` in one snapshot was missing `ormeta`/`ormidx`
registrations that existed in the actually-running package). Always ask
for current versions when Fortran-level behavior is in question, rather
than reasoning from a possibly-outdated copy.

## 1. Architecture

```
orm.fit()
  |-- ordinary (non-clustered) fixed-effects warm-up via ormfit()
  |     (also used standalone for the intercept-only null deviance,
  |      the offset-only deviance, etc.)
  |
  |-- if cluster present and p>0: refit a clustered, covariate-free
  |     NULL model too (needed for a correct LR test -- see Sec. 5),
  |     itself warm-started from a fresh, non-clustered intercept-only
  |     refit (see Sec. 7's "bad warm start" bug)
  |
  '-- ormrfit()  [the actual random-effects fitting entry point]
        |
        '-- outer loop (SQUAREM acceleration, Varadhan & Roland 2008
              "SqS3" scheme), each cycle calling onestep() 2-3 times
              |
              '-- onestep(), per call:
                    1. clusterModeFind(): per-cluster mode (gamma_c)
                       and curvature (kappa_c) of the random effect,
                       given current (alpha, beta, sigma/sigma1/sigma2)
                    2. agqStep(): AGQ-integrated log-likelihood,
                       gradient, Hessian at the current parameters
                    3. Newton step:
                       - mre INACTIVE: joint step for (alpha, beta)
                         only; sigma gets its OWN separate 1-D profile
                         search via optimize() (see Sec. 4)
                       - mre ACTIVE: joint step for
                         (alpha, beta, log(sigma1), sigma2) together
                    4. step-halving surrogate (EM/MM-style trial
                       objective, cheap ormeta-only evaluation)
              after outer-loop convergence: an OPTIONAL nAGQ escalation
              check (see Sec. 6 for a crash this exposed)
        |
        '-- sparseInfoMatrix() / sparseMissingInfo(): final, corrected
              (Louis's-identity-adjusted) information matrix
```

`infoMxop()` operates on the assembled information matrix afterward
(inversion, extracting `sigma1`/`sigma2` covariance, etc.) — see
Sec. 4 of the methodology doc's testing discussion and Sec. 3 below for
how it locates the random-effect-scale elements.

Data structures `ia`/`ia2`/`sgn`/`ib`/`nb` are a pure function of
`y`/`y2`/`k` (computed once via `ormidx`, reused across every call).
`theta_full`'s layout differs by whether `mre` is active:
`mre` inactive: `c(alpha, beta, log(sigma))`;
`mre` active: `c(alpha, beta, log(sigma1), sigma2)` (no separate `tau`
at all in the active case).

## 2. Rejected designs for the weighted random effect

Three designs were tried, in order, before settling on the additive
`(sigma1, sigma2)` form described in the methodology doc.

**Design 1 — `(wl, rho)`**: `w(tre) = wl + (1-wl)*rho^tre`, two
estimated parameters. Abandoned: genuinely unstable through direct
testing — an early Newton step could drive `rho` to the exact
floating-point boundary, collapsing `wl`'s own curvature. A grid-scan
of the likelihood confirmed a genuine competing-mode/weak-identifiability
problem, not a fixable numerical issue. Curvature at the true values
scaled with `nc` (did not plateau), but the multiple-mode risk was
real and persistent regardless of `nc`.

**Design 2 — `(sigma, re_mult)`**: `(1 - mre*re_mult)*u_i`,
`u_i ~ N(0, sigma^2)`. Linear and well-identified *in isolation*
(confirmed via direct curvature check). Failed in practice: when `mre`
is active and `re_mult>1`, `sigma -> 0` while `re_mult -> -infinity`
*together* — confirmed via a direct log-likelihood grid scan showing
the drift was in the *wrong* direction (the true parameter values gave
a clearly higher likelihood than the fitted, drifted ones). Root cause:
`re_mult`'s scale is defined only *relative to* `sigma`, so as
`sigma -> 0`, `re_mult` becomes structurally meaningless. Folding
`tau = log(sigma)` into the same joint Newton step as `re_mult` (rather
than a separate profile search) reduced but did not eliminate the
co-drift.

**Design 3 — `(sigma1, sigma2)` [current, chosen]**: each an *absolute*
scale, not a ratio relative to the other. Curvature check: correlation
between `sigma1` and `sigma2` stayed in the range -0.05 to -0.12 across
every scenario tested (vs. 0.7-0.8 for Design 2's `(sigma, re_mult)`),
with both eigenvalues strongly positive even at small `nc` (e.g. 22).
This is *why* the additive form works where the multiplicative one
didn't: each parameter is identified by its own share of the data
variance, with no ratio relationship to fail as either one moves.

`sigma2`'s starting value must be exactly `0` (not `sigma.init`, not
any other nonzero value) — confirmed critical. Starting at a nonzero
value biased the fit toward that sign persistently, independent of
sample size, since `sigma2`'s own gradient
(`d_sigma2 = mre * gamma_node_c[cluster]`, no chain-rule scale factor)
is well-defined and correctly signed at exactly zero, letting the
Newton step discover the correct sign directly from the data.

## 3. `infoMxop`: dropping positional assumptions

Originally, `infoMxop` located the random-effect scale element(s) by
assuming `xname`'s *last* entry was named `'log(sigma)'`
(`nsigma <- length(xname) && xname[length(xname)]=='log(sigma)'`).
This was generalized to *search* `xname` for `'log(sigma)'` and
`'sigma2'` by name, wherever they occur, with `nsigma` becoming a plain
count (0, 1, or 2) rather than a boolean:

- `i='log(sigma)'` renamed to `i='sigma_parameters'`; returns a
  scalar variance when only `log(sigma)` is present, or a 2x2
  covariance matrix (with dimnames `c('log(sigma)','sigma2')`, in that
  fixed order regardless of `xname`'s actual order) when both are
  present.
- The boundary-failure NA handling generalized: a scalar `NA` or an
  all-`NA` 2x2 matrix, matching whether one or two scale parameters
  are present, with an explanatory warning either way.
- The ridge-retry fallback (for a genuinely singular info matrix at a
  `sigma`-near-zero boundary) generalized from ridging a single assumed
  position (`nv+1`) to ridging *all* found sigma-related positions.
- `orm.rfit.r`'s own `sparseInfoMatrix` was changed to emit
  `'log(sigma)'` (not `'log(sigma1)'`) as the label even when `mre` is
  active, so `infoMxop`'s search finds it consistently regardless of
  whether `mre` is in use.

All of the above was verified against an independent, direct matrix
inversion (not just internal consistency), including a check with
`sigma2` appearing *before* `log(sigma)` in `xname`, to confirm the
search is genuinely order-independent.

## 4. The joint-Newton-step regression and revert

Folding `sigma` into the *same* joint Newton step as `(alpha, beta)`
(rather than giving it a separate 1-D profile search) was originally
added to fix Design 2/3's negative-`sigma2`-sign-bias problem (Sec. 2).
It was then generalized to run for *all* cases, including `mre`
inactive — **this was a real regression**, caught only via a real
user dataset (`nc=250` clusters of size 2, continuous-ish `y`): the fit
never converged, with `sigma` drifting slowly and monotonically for
100 outer iterations without stabilizing (originally converged in 3-24
outer iterations before the fold).

**Cause**: a single linearized Newton step for `log(sigma)` can
systematically undershoot when `sigma`'s own curvature is asymmetric —
the original, separate `optimize()`-based profile search doesn't have
this problem because it's a direct, non-linearized 1-D search at each
outer iteration.

**Fix**: reverted to the separate profile search specifically when
`mre` is inactive (matching the original, pre-fold design); kept the
joint step only for the `mre`-active case, where it was actually
validated to fix a genuine, different problem. Verified: the original
`nc=250` regression now converges in 12 outer iterations; the
`mre`-active negative-`sigma2` scenario that motivated the joint step
in the first place still converges correctly.

**Lesson**: a fix validated for one code path (here, `mre` active)
should not be generalized to *all* paths without re-validating each
one — the two cases have genuinely different curvature behavior for
`sigma`, and treating them identically broke the far more common case.

## 5. The likelihood-ratio-test bug (`orm.fit.r`)

`orm.fit` accumulates a `loglik` vector as it runs through its fitting
stages (intercepts-only, intercepts+offset, intercepts+x, ...), and the
stats block computes `model.lr <- loglik[length(loglik)-1] -
loglik[length(loglik)]`, i.e. it assumes the *last two* elements are
`[null, full]`.

**Bug**: when `cluster` was added, only the *full* clustered deviance
was appended (`loglik <- c(loglik, zr$loglik)`) — no clustered null was
ever computed. So `loglik[length(loglik)-1]` silently grabbed the
**non-clustered full-model deviance** instead of any null at all. The
resulting "Model L.R." was actually testing "does the random effect
improve fit over the fixed-effects-only model" — a test related to
`sigma`, not to `beta` — which is why it bore no relationship to
`beta`'s own Wald chi-square (confirmed: Wald `chi-sq ~ 1024`
vs. reported LR `chi-sq` of `177`/`224` for the two models).

**Fix**: when `cluster` is present and `p>0`, also fit a clustered
*null* model (intercepts + the same random-effects structure, no `x`)
and append **both** its deviance and the full model's deviance. This
keeps the existing "last two elements are `[null, full]`" indexing
correct without needing to touch the stats block itself. Verified
directly: the corrected Model L.R. (`121.6`) now agrees closely with
the Wald chi-square (`115.7`) on a real test case, rather than
differing by 5-9x.

**Side effect worth knowing**: every element pushed into `loglik` is
now also *named* at the point it's computed
(`'intercepts'`, `'intercepts+offset'`, `'intercepts+x'`,
`'intercepts+random effects'`, `'intercepts+x+random effects'`),
regardless of which combination of censoring/offset/covariates/
clustering branches actually ran — see the `deviance` element's own
`@returns` documentation in `orm.fit.r` for the full, ordered
composition.

## 6. Bug log: numerical robustness (symptom -> cause -> fix)

Format: what was observed -> what was actually wrong -> what changed.
Kept terse and lookup-oriented; see git history / commit comments in
the actual source for full narrative detail on any of these.

**`clusterModeFind failed (ormeta salloc=999)` under `probit`, first
outer iteration, `sigma=1`.**
Cause: `ormeta.f90`'s own documentation states `salloc=999` "signals
the caller to step-halve, same as `ormll`'s `what=1` path" — a
*designed* retry signal, not a hard failure. `ormll`'s own caller (the
ordinary NR loop) already honors this; `clusterModeFind` did not — any
nonzero `salloc` was treated as immediate, unrecoverable failure.
Fix: proper step-halving retry loop in `clusterModeFind`, evaluating
`ormeta` once at the current point, then halving the Newton step (down
to a small minimum) on `salloc != 0` before giving up. Caveat: if the
*very first* evaluation (no random effect contribution yet) already
fails, there's no prior good point to halve back toward — that case
still fails immediately, by design, since step-halving cannot rescue an
already-invalid starting point.

**Why `probit` specifically triggers `salloc=999` more than `logit`.**
`ormll.f90`'s density formulas: logistic's `pdf = f*(1-f)` is computed
*from the CDF value* and decays gracefully; probit's `pdf =
exp(-x^2/2)/sqrt(2*pi)` is a direct Gaussian tail that underflows to
exactly `0.0` at a much smaller `|x|`. Two adjacent intercepts
perfectly distinguishable (`d>0`) under logistic can both round to the
same saturated CDF value under probit, making their difference (a
category probability) exactly zero.

**`Error in eval_ormeta(gamma_c_new) : NA/NaN/Inf in foreign function
call (arg 5)`.**
Cause: R's own `.Fortran()` interface performs a hard, unconditional
check for `NA`/`NaN`/`Inf` in numeric arguments *before* dispatching to
Fortran at all, and throws an **uncatchable** R-level error if it finds
one — before `salloc` is ever produced, so the step-halving fix above
never even got a chance to run. `gamma_c_new` itself was going
non-finite (a large, still-growing `sigma` driving `curv` toward its
`-1e-8` floor, producing an enormous, eventually infinite, Newton
step).
Fix: explicit `is.finite()` checks immediately before *every*
`.Fortran(F_ormeta,...)` / `.Fortran(F_ormll,...)` call in
`orm.rfit.r` (`clusterModeFind`, both of `agqStep`'s loops,
`sparseMissingInfo`, the step-halving surrogate). A non-finite trial
point is treated exactly like `salloc=999` (or, in the surrogate,
returns `Inf` to match its existing graceful-rejection convention),
rather than reaching `.Fortran()` and crashing.

**cloglog's asymmetric fragility (the `nAGQ` escalation check
specifically).**
`ormll.f90`'s `cloglog` density (`case(4): p = exp(x - exp(x))`) is a
one-sided, double-exponential decay — measured directly: it underflows
to exact `0.0` by `x ~ 7`, versus `x ~ 37` for logistic and `x ~ 39`
for probit (over 5x smaller threshold). The lower tail behaves like
logistic's (no fragility there) — this is specifically an upper-tail,
asymmetric phenomenon. Gauss-Hermite nodes spread wider as `nAGQ`
increases; escalating from `nAGQ=7` to `11` pushes the widest nodes far
enough that, combined with `sigma` scaling, some observation's shifted
linear predictor can cross that much-smaller cloglog threshold. This
is a genuine, structural property of the link, not a bug — the fit
itself converges fine at the smaller `nAGQ`; only the *optional*
wider-`nAGQ` stability check is affected, and it now fails gracefully
(see next entry) rather than crashing.

**`Error in if (rel.change < nAGQ.tol) break : argument is of length
zero`.**
Cause: the `nAGQ` escalation check (run *after* the outer loop already
converged, to verify a larger `nAGQ` wouldn't change the answer) never
validated whether its own `clusterModeFind`/`agqStep` calls succeeded.
`agqStep`'s failure return, `list(fail=TRUE, code=...)`, has no
`$loglik` element at all, so `$loglik` on it is `NULL`; `abs(NULL -
ll.cur)` silently collapses to `numeric(0)`, and
`if(numeric(0) < nAGQ.tol)` throws exactly this error, with no
indication of the real cause. This is the failure mode the cloglog
entry above actually hits.
Fix: explicit `$fail` check after each of the three calls
(`clusterModeFind`, `agqStep` at `cur.nAGQ`, `agqStep` at `next.nAGQ`).
Since this whole block is optional verification (the fit already
converged before reaching it), any failure is treated as "cannot
verify a larger `nAGQ` would help" and the loop simply breaks, keeping
the already-converged result at the current `nAGQ`.

**Bad warm start for the clustered null-model refit (added for the
LR-test fix in Sec. 5).**
Symptom: under `probit`, the null-model refit itself failed (zero NR
iterations logged — failing at the very starting point, not from an
overshooting step). A second attempt, warm-starting from `kof[1:k]`
(the full model's own intercepts) still failed the same way, just
moved to a different call.
Cause: `kof[1:k]` are the full model's intercepts, spread out to
*compensate* for `p` covariates. Removing the covariates while keeping
those same, compensated intercepts is often invalid: measured directly
in a synthetic test (large covariate effect, `k=30` intercepts) — the
smallest adjacent-category probability gap shrank from `0.0025`
(properly-refit intercepts) to `0.00037` (naive, covariate-compensated
intercepts) — nearly 7x smaller, and squarely in the range where
probit/cloglog underflow to an exact zero.
Fix: use `finverse(pp)` — the same link-appropriate, purely
marginal-frequency-based intercept-only starting guess `orm.fit`
already computes for its own very first, non-clustered fit — recomputed
fresh at the point of the null-model refit (not reused from whatever
`initial` happened to be, so it's robust even with a user-supplied
custom `initial`).

## 7. Debug infrastructure: `Fdebug`

`Hmisc::Fdebug(opt)` returns either a real print-with-label function or
a do-nothing function, decided once at creation time based on
`getOption(opt, FALSE)`:

```r
deb <- Hmisc::Fdebug('orm.re.debug')
...
deb(gamma_c)   # prints only if options(orm.re.debug=TRUE) is set;
               # otherwise a true no-op -- the argument is never
               # forced, so this costs almost nothing even for an
               # expensive expression
```

Labels come from `deparse(substitute(x))` (the object's own expression
text) plus the calling function's name (via `sys.call(-1)`, captured at
`Fdebug()` creation time — correctly resolves to the function `deb` was
created inside, even when `deb()` is later invoked from a nested
closure like `eval_ormeta`).

Currently wired into `clusterModeFind` (on `gc`/`lp`/`w$salloc`/`alpha`)
and at two checkpoints in `ormrfit`'s and `orm.fit`'s own top-level
setup (`p`/`dim(x)`/`kof`/`beta`/`xbeta`/`base_lp`), added specifically
to trace the warm-start bug in Sec. 6.

**Explicit correction to keep in mind**: an earlier version of this
mechanism also tied `ormeta`'s own internal (Fortran-level, `intpr()`)
debug flag to the same option (`ormeta_debug <- as.integer(...)`,
passed as `debug=ormeta_debug`). This was deliberately removed at the
developer's request — `ormeta`'s own `debug` argument should stay a
plain, hardcoded `0L`; only `deb()` calls should be gated by the
option. Keep these two concerns separate if extending this mechanism
further.

## 8. Validation methodology used throughout

Worth preserving as a template for future changes to this code:

- **Brute-force numerical comparison** for the base single-`sigma`
  case (an independent, from-scratch computation to check against,
  not just internal self-consistency).
- **Direct curvature/correlation checks** for identifiability claims
  (Sec. 2) — computed the actual Hessian/correlation numbers rather
  than reasoning about them abstractly.
- **Synthetic simulations at varying scale** — the joint-step
  regression (Sec. 4) was only caught by testing at `nc=250`,
  cluster-size 2 (a scale earlier test scenarios, with small numbers
  of larger clusters, never exercised).
- **Mocked `.Fortran` stand-ins** for control-flow-only tests when the
  real compiled routine isn't available in-session (e.g. testing
  `clusterModeFind`'s step-halving retry logic against a synthetic
  `ormeta` that deliberately fails in a controlled way).
- **Reading the actual Fortran source, not guessing from symptoms
  alone**, once available — several diagnoses (the `salloc=999`
  contract, the exact `pdf`/`cdf` formulas, cloglog's asymmetric decay
  rate) came from reading `ormll.f90`/`ormeta.f90` directly rather than
  inferring behavior from R-level symptoms.

## 9. Known limitations / open items for future work

**Y-dependent effects (yde) / partial proportional odds — needs
dedicated examination before combining with `cluster()`/`mre`.** See
the methodology doc's Sec. 6 for the statistical motivation; the
specific technical questions for a future session:

1. `ormeta`'s signature already carries `lp1`/`lp2` (potentially
   different linear predictors depending on which side of a category
   cutpoint an observation's `y` falls on) and separate `ia`/`ia2`
   index arrays for exactly this purpose. **Every current call site in
   `orm.rfit.r` sets `lp1=lp2=lp` unconditionally** (confirmed by
   direct grep across the file) — full proportional odds is baked in
   everywhere the random-effects code touches `ormeta`/`ormll`.
2. `clusterModeFind`'s mode-finding score/curvature formulas
   (`Sg <- rowsum(wtre*w$g, cluster,...)`, `Sh <- rowsum(wtre^2*w$h,
   cluster,...)`) currently sum a single `g`/`h` per observation. If a
   y-dependent covariate makes an observation's effective linear
   predictor differ depending on cutpoint, does a single scalar `g`/`h`
   per observation still correctly represent that observation's
   contribution to the cluster-level mode-finding score/curvature, or
   does the derivation need to be redone from the partial-PO
   likelihood directly?
3. `agqStep`'s pseudo-covariate/Hessian-augmentation logic for
   `sigma1`/`log(sigma)`/`sigma2` (the `d_logsigma1`/`d_sigma2` chain-rule
   terms) was derived under the assumption of a single, shared linear
   predictor per observation. Does this generalize cleanly, or does a
   y-dependent model need its own analogous derivation?
4. `sparseMissingInfo`'s missing-information (Louis's identity)
   correction is built from `ia`/`ia2` under the assumption that these
   index a single set of per-observation quantities. Y-dependent
   effects would make some of those quantities themselves
   cutpoint-dependent — does the current sparse correction still apply
   as-is, or does it need extension?

**Other open items, not yet addressed:**

- Partial identifiability of `(sigma1, sigma2)` at small `nc` — noted
  early in development as an open design question (a
  sample-information-based fallback to a `sigma`-only model was
  discussed but never implemented); the `mre`-constant-variation check
  handles only the fully-degenerate case, not partial identifiability
  at small sample sizes.
- `lpe` (per-observation likelihood contribution) is not computed for
  clustered fits — currently returns `NA`.
- Fitting speed / SQUAREM cycle cost as `p` or the number of
  intercepts `k` grows very large has not been specifically profiled.
- The `nAGQ`-escalation fragility under `cloglog` (Sec. 6) is
  documented and handled gracefully, but not itself "fixed" — a real
  fix, if ever wanted, would need an asymmetric AGQ node grid (denser
  toward the link's fragile tail) rather than the current symmetric
  one.

## References

- Varadhan, R. and Roland, C. (2008). Simple and globally convergent
  methods for accelerating the convergence of any EM algorithm.
  *Scandinavian Journal of Statistics* 35(2), 335-353. — the "SqS3"
  SQUAREM scheme used for the outer-loop acceleration.
