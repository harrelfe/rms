#' Ordinal Regression Model Fitter
#'
#' @description
#' Fits ordinal cumulative probability models for continuous or ordinal
#' response variables, efficiently allowing for a large number of
#' intercepts by capitalizing on the information matrix being sparse.
#' Five different distribution functions are implemented, with the
#' default being the logistic (yielding the proportional odds model).
#' Penalized estimation and weights are also implemented, as in
#' [lrm.fit()]. The optimization method is Newton-Raphson with
#' step-halving, or the Levenberg-Marquardt method. The latter has
#' been shown to converge better when there are large offsets, or when
#' the likelihood surface is moderately ill-conditioned (see
#' **Details**). Execution time is fast even for hundreds of thousands
#' of intercepts. The limiting factor is the number of intercepts
#' times the number of columns of `x`.
#'
#' `orm.fit` also fits a single-level random-intercept extension of
#' this model when `cluster` is specified, for clustered or repeated-
#' measures data, optionally with a user-specified, continuously- or
#' discretely-varying weight on the random effect (see `mre` below and
#' **Details**). See **Details** for a full description of the
#' random-intercept model, how it is fit, and how sparseness is
#' preserved even with thousands of intercepts.
#'
#' @details
#' ## Random-intercept model
#' When `cluster` is given, `orm.fit` fits `k` ordinal intercepts and
#' `p` regression coefficients as usual, plus a normally-distributed
#' random intercept `gamma_j ~ N(0, sigma^2)` for each of the distinct
#' clusters (levels of `cluster`). Fitting integrates out the random
#' effects to obtain the marginal likelihood, using adaptive
#' Gauss-Hermite quadrature (AGQ) rather than a Laplace approximation.
#' This choice was empirically driven: extensive validation against
#' `ordinal::clmm2` showed that Laplace approximation (equivalent to
#' `nAGQ=1`) produces substantial, systematic bias in both the fixed
#' effects and, especially, the variance component -- worsening as
#' cluster informativeness decreases (small clusters, large `sigma`),
#' which is the regime most common in practice. `nAGQ` is not held
#' fixed either: it is escalated automatically, refitting (warm-
#' started) at each successively larger value in `nAGQ.grid` until the
#' log-likelihood stabilizes within `nAGQ.tol`, since no single fixed
#' value of `nAGQ` was found to be adequate across the full range of
#' `sigma` a user might encounter.
#'
#' ## Weighted random effects via `mre`
#' When a fixed, user-specified per-observation vector `mre` is also
#' given, the random intercept's own contribution to the linear
#' predictor is no longer a constant `gamma_j` -- it is scaled by a
#' quantity that can vary, continuously if desired, across an
#' individual cluster's own observations:
#' \deqn{\eta_{ij} = \ldots + [\sigma_1(1-\mbox{mre}_{ij}) +
#'   \sigma_2\,\mbox{mre}_{ij}]\, v_i}
#' where `v_i ~ N(0,1)` is a standardized, per-cluster random effect
#' and `sigma1`, `sigma2` are two estimated scale parameters --
#' `sigma1` constrained positive via a `log` link (to pin down `v_i`'s
#' otherwise arbitrary sign, since flipping the sign of both `sigma2`
#' and every `v_i` simultaneously leaves the likelihood unchanged
#' unless something anchors it), `sigma2` left unconstrained, so it
#' may legitimately come out negative. `mre` is entirely user-computed
#' and fixed -- it is data, exactly like a column of `x`, never
#' estimated -- while `sigma1` and `sigma2` are estimated jointly with
#' `alpha` and `beta` in the same Newton step. `mre=NULL` (the
#' default) recovers the ordinary, single-`sigma` model exactly, with
#' no loss of speed or accuracy.
#'
#' This additive form -- `sigma1` and `sigma2` each multiplying a
#' separately-known quantity, `(1-mre)v_i` and `mre \, v_i`
#' respectively -- was chosen over an earlier, multiplicative design
#' (a single `sigma` combined with a separately estimated multiplier)
#' specifically because it keeps `sigma1` and `sigma2` close to
#' orthogonal. Extensive testing found the multiplicative form's extra
#' parameter became structurally unidentifiable whenever the base
#' `sigma` approached its lower boundary of zero, causing the fit to
#' converge to an arbitrary, meaningless value for it; `sigma1` and
#' `sigma2` remain well-behaved instead (their estimated correlation
#' stayed well under `0.2` in magnitude across every scenario tested,
#' including deliberately small samples and scale parameters of
#' opposite sign), since each is an absolute scale directly identified
#' by its own share of the data rather than a ratio relative to the
#' other.
#'
#' `mre` must vary *within* at least some clusters for `sigma1` and
#' `sigma2` to be separately identifiable; `orm.fit` checks this
#' directly and stops with an informative error if `mre` has no
#' variation anywhere in the data (e.g. is constant at `0`, at `1`, or
#' at any other single value for every observation). A cluster
#' contributing only a single observation is harmless and needs no
#' special handling -- it simply contributes nothing toward separating
#' `sigma1` from `sigma2`, the same way a covariate with no variation
#' in one subgroup does not hurt an ordinary regression, as long as
#' *some* clusters elsewhere in the data do have within-cluster `mre`
#' variation.
#'
#' A natural and recommended choice for a Markov-1 model on repeated
#' ordinal responses (a proportional-odds model whose linear predictor
#' includes a term for the previous response) is a piecewise-constant
#' `mre`: `0` at each subject's first follow-up time and `1` at every
#' subsequent one (see **Examples**). The motivation is that a plain,
#' constantly-weighted random effect double-counts its own influence
#' once a lag term is present -- the previous response already carries
#' the random effect's influence forward, so adding it again at full
#' strength at every subsequent visit can induce an unrealistically
#' large or ever-growing correlation across time. Letting `sigma2`
#' differ from `sigma1` -- including, if the data call for it, being
#' of the opposite sign -- lets the random effect's influence after
#' the first visit either add to or subtract from the correlation the
#' lag term alone would induce, as the data determine. `mre` is not
#' restricted to this `0`/`1` step function, however; any fixed,
#' user-computed function of the data (continuously varying with
#' follow-up time, for example) is equally valid.
#'
#' ## Likelihood ratio tests under clustering
#' When `cluster` is given together with `p>0` covariates, `orm.fit`
#' fits an *additional*, otherwise-unnecessary clustered model:
#' intercepts plus the same random-effects structure (including `mre`,
#' if given), but without `x`. This provides a clustered null deviance
#' that is directly comparable to the clustered full model's deviance,
#' so that the reported model likelihood ratio chi-square is a genuine
#' test of the covariates given the random-effects structure already
#' in the model -- not a test confounded with whether the random
#' effect itself improves the fit. See the `deviance` element of
#' **Returns** for the full, named composition of the returned
#' deviance vector.
#'
#' ## Inner and outer loops
#' Fitting uses an outer, block-alternating loop, each iteration of
#' which does the following:
#'
#' 1. **Per-cluster mode-finding** (the inner loop): a scalar
#'    Newton-Raphson iteration finds each cluster's own posterior mode
#'    and curvature for its random intercept, at the current
#'    `(alpha, beta, sigma)`. This is vectorized across *all* clusters
#'    simultaneously -- one lean Fortran call per Newton step, not one
#'    call per cluster -- so its cost is `O(n)` regardless of the
#'    number of clusters.
#' 2. **AGQ-reweighted Newton step for `(alpha, beta)`**: at each of
#'    the current `nAGQ` quadrature nodes (centered and scaled per
#'    cluster using the mode and curvature just found), the same,
#'    completely unmodified `orm.fit` likelihood/gradient/Hessian
#'    machinery is called once, with that node's random-intercept
#'    value folded into the offset and its posterior node weight
#'    folded into the case weights. The `nAGQ` results are summed to
#'    form one Newton step. Because this reuses `orm.fit`'s own
#'    tri-band (or, under interval censoring, general sparse triplet)
#'    intercept-Hessian representation and its `O(k)` solve unchanged,
#'    the random-intercept extension inherits the same scalability to
#'    very large numbers of intercepts that `orm.fit` has without
#'    random effects -- it simply pays this cost `nAGQ` times per
#'    outer iteration instead of once.
#' 3. **Fitting the random-effect scale parameter(s)**. When `mre` is
#'    not given, this is a one-dimensional profile search for `sigma`
#'    via `optimize()`, using only the lean, `O(n)` per-observation
#'    routine from step 1 (no gradient or Hessian needed) -- cheap
#'    relative to step 2. This separate-search design was tried for
#'    the `sigma1`/`sigma2` case too, but found to converge to a
#'    persistent wrong value whenever the true `sigma2` had the
#'    opposite sign from its (necessarily sign-unknown-in-advance)
#'    starting value, regardless of sample size -- confirmed to be a
#'    genuine optimization-path problem, not a sign of poor
#'    identifiability, since a direct grid-scan of the likelihood
#'    always found a well-behaved optimum near the true values. When
#'    `mre` *is* given, `sigma1` and `sigma2` are therefore instead
#'    estimated jointly with `alpha` and `beta` in the very same
#'    Newton step as step 2, `sigma2` started at exactly `0` (not an
#'    arbitrary nonzero value) so the step itself discovers the
#'    correct sign directly from the data rather than being biased
#'    toward whichever sign it happened to start on.
#'
#' Because this scheme converges only linearly near the optimum
#' (unlike `orm.fit`'s own Newton-Raphson, which converges
#' quadratically), its convergence criteria differ correspondingly:
#' the gradient tolerance is used as an absolute, not an `n`-scaled,
#' criterion (an `n`-scaled tolerance is nearly a no-op at realistic
#' sample sizes for a scheme with only linear convergence), and a
#' separate allowance accepts convergence once `sigma`'s absolute
#' change falls below a small floor even if its *relative* change has
#' not, since relative changes become numerically meaningless once
#' `sigma` is already effectively at its lower boundary of zero.
#'
#' The outer loop's own convergence is accelerated using SQUAREM
#' (Varadhan and Roland, 2008), a generic extrapolation method for
#' fixed-point iterations that share this scheme's linear/geometric
#' convergence pattern. Every SQUAREM cycle runs two ordinary outer
#' iterations, uses them to estimate the local geometric contraction
#' rate, and extrapolates a candidate point that would take many more
#' ordinary iterations to reach directly. That candidate is never
#' accepted on faith: it is itself run through one more ordinary outer
#' iteration ("stabilization"), and is kept only if that step's own
#' parameter change is no larger than the plain, unaccelerated
#' iteration's would have been -- so acceleration can only match or
#' improve on plain iteration, never meaningfully underperform it. The
#' extrapolation's step length is additionally bounded by an adaptive
#' cap that starts conservative and grows after each accepted
#' extrapolation, shrinking back after a rejected one -- an unbounded
#' step length is only trustworthy once the local trajectory has
#' settled into a clean geometric pattern, and can otherwise produce
#' wild, wastefully expensive trial points. In testing this combination
#' reduced wall clock fitting time by roughly 15-20% relative to
#' unbounded extrapolation, and reduced the number of outer iterations
#' needed by roughly 3-to-7-fold relative to no acceleration at all, on
#' both a very small (20-observation, 10-cluster) and a more
#' typical-scale (400-observation, 200-cluster) dataset, while landing
#' at coefficient estimates matching plain (unaccelerated) iteration to
#' well within the convergence tolerance itself.
#'
#' ## Preserving sparseness
#' Every piece of the random-intercept computation was designed so
#' that its cost scales with `n`, the number of clusters, or cluster
#' size -- never with `k^2` or `k^3`:
#'
#' - The per-cluster mode-finding and AGQ node evaluation touch `k`
#'   only as a lookup into `alpha`, never as a matrix dimension.
#' - The `(alpha, beta)` Newton step reuses `orm.fit`'s existing
#'   tri-band (no interval censoring) or general sparse triplet
#'   (interval censoring present) representation of the intercept
#'   Hessian and its `O(k)` tridiagonal solve, completely unchanged.
#' - The final, corrected information matrix's missing-information
#'   correction (via Louis's (1982) identity, which corrects for the
#'   variance a fixed-posterior-weight Hessian omits) is accumulated
#'   as a sparse triplet list: each cluster contributes only to the
#'   alpha indices *that cluster's own observations actually touch* --
#'   at most about twice its own size, however large `k` is. The
#'   total number of nonzero entries this produces is bounded by
#'   roughly the number of clusters times the square of a typical
#'   cluster's size, independent of `k`, for realistic (bounded
#'   cluster size) designs -- confirmed empirically (e.g. 588 nonzero
#'   entries out of 11,175 possible for a `k=149`, 75-cluster
#'   example). This lets `orm.fit` fit random-intercept models with
#'   thousands of intercepts in about a second, versus tens of seconds
#'   to minutes for `ordinal::clmm2` on the same data.
#' - `log(sigma)`'s standard error (or the covariance of `log(sigma)`
#'   and `sigma2`, when `mre` is given) is obtained by extending only
#'   the non-intercept (`b`/`ab`) pieces of `info.matrix` by one or two
#'   rows and columns -- the intercept-intercept block is never
#'   touched by this extension, so [infoMxop()] and other code that
#'   consumes `info.matrix` require no special handling for the
#'   random-effects case beyond what interval censoring already
#'   required. [infoMxop()] locates these elements by searching
#'   `info.matrix$xname` for `"log(sigma)"` and `"sigma2"` wherever
#'   they occur, rather than assuming a fixed position, so their
#'   presence and order elsewhere in `xname` does not matter.
#'
#' ## Other computational details
#' - `opt_method='LM'` uses an *asymmetric* Levenberg-Marquardt damping
#'   schedule (`lambda` divided by 1.5 on an accepted step, multiplied
#'   by 10 on a rejected one) rather than a symmetric `/10`/`*10` rule.
#'   The symmetric rule can enter a permanent limit cycle on
#'   moderately ill-conditioned (but not literally singular)
#'   likelihood surfaces -- alternating forever between a damping
#'   level that barely succeeds and one that barely fails, with the
#'   true optimum never actually reached. The asymmetric rule was
#'   confirmed, by direct inspection of `lambda`'s trajectory on such
#'   a case, to break this cycle, at no cost in iterations for
#'   well-conditioned problems.
#' - Fitted standard errors (via [infoMxop()] operating on
#'   `info.matrix`) always come from this corrected, not naive,
#'   information matrix when random effects are present -- ignoring
#'   the missing-information correction was found in validation
#'   studies to understate intercept standard errors by as much as
#'   40%.
#' - Requesting `i='sigma_parameters'` from [infoMxop()] can encounter
#'   a genuinely singular information matrix when a fit has converged
#'   with `sigma` (or `sigma1`) at or very near its lower boundary of
#'   zero -- at that boundary its information content correctly
#'   collapses to zero (a known, non-regular feature of variance-
#'   component estimation, not a numerical defect), and `infoMxop`
#'   returns `NA` (a scalar, or an all-`NA` `2x2` matrix when `sigma2`
#'   is also present) with an explanatory warning rather than raising
#'   a raw singularity error in that case.
#' - The effective-sample-size adjustment for censored observations
#'   (contributing to the `"ESS"` element of `stats`) is not yet
#'   available when `cluster` is specified; `stats["ESS"]` reflects
#'   only the uncensored-observation count in that case, and a warning
#'   is issued.
#'
#' @param x design matrix with no column for an intercept
#' @param y response vector, numeric, factor, or character. The
#'   ordering of levels is assumed from `factor(y)`.
#' @param cluster optional vector, the same length as `y`, of cluster
#'   identifiers (e.g. subject ID for repeated measurements). When
#'   given, `orm.fit` fits an additional normally-distributed random
#'   intercept for each distinct cluster; see **Details**. Left as
#'   `NULL` (the default) to fit the ordinary fixed-effects model.
#' @param mre optional fixed, numeric vector the same length as `y`,
#'   used only when `cluster` is also given: a user-computed
#'   "multiplier for random effects" governing how the random
#'   intercept's own contribution is weighted at each observation, via
#'   two estimated scale parameters `sigma1` and `sigma2` in place of
#'   the ordinary single `sigma` (see **Details**). Left as `NULL`
#'   (the default) to fit the ordinary, single-`sigma` random-
#'   intercept model. For a Markov-1 model on repeated ordinal
#'   responses, a natural choice is `mre` equal to `0` at each
#'   subject's first follow-up and `1` at every later one (see
#'   **Examples**), but any fixed function of the data, continuous or
#'   discrete, is valid. `mre` must vary within at least some
#'   clusters; a completely constant `mre` (whether `0`, `1`, or any
#'   other single value across the whole dataset) is rejected with an
#'   informative error, since it would leave `sigma1` and `sigma2`
#'   unidentifiable.
#' @param family a character value specifying the distribution family,
#'   corresponding to logistic (the default), Gaussian, Cauchy, Gumbel
#'   maximum (`exp(-exp(-x))`; extreme value type I), and Gumbel
#'   minimum (`1-exp(-exp(x))`) distributions. These are the
#'   cumulative distribution functions assumed for
#'   `Prob[Y >= y | X]`. The `family` argument can be an unquoted or a
#'   quoted string, e.g. `family=loglog` or `family="loglog"`. To use
#'   a built-in family, the string must be one of the following,
#'   corresponding to the previous list: `logistic`, `probit`,
#'   `loglog`, `cloglog`, `cauchit`.
#' @param offset optional numeric vector containing an offset on the
#'   logit scale
#' @param initial vector of initial parameter estimates, beginning
#'   with the intercepts. If `initial` is not specified, the function
#'   computes the overall score \eqn{\chi^2} test for the global null
#'   hypothesis of no regression. `initial` is padded to the right
#'   with zeros for the regression coefficients, if needed. When
#'   censoring is present, `initial` can also be a list with elements
#'   `time` and `surv` from the `npsurv` attribute of the `y` element
#'   of a previous fit. This is useful when bootstrapping, for
#'   example. When `cluster` is given, `orm.fit`'s own fixed-effects
#'   initial-value refinement (run once, regardless of `cluster`) is
#'   reused as the warm start for the random-intercept fit, so
#'   `initial` need not account for random effects.
#' @param opt_method set to `"LM"` to use Levenberg-Marquardt instead
#'   of the default Newton-Raphson; see **Details** for the asymmetric
#'   damping schedule used.
#' @param maxit maximum no. iterations (default=`30`). For the
#'   random-intercept outer loop this is instead governed by
#'   `maxit.outer`.
#' @param eps difference in \eqn{-2 \log} likelihood for declaring
#'   convergence. Default is `.0005`. This handles the case where the
#'   initial estimates are MLEs, to prevent endless step-halving.
#' @param gradtol maximum absolute gradient before convergence can be
#'   declared. For the fixed-effects fit, `gradtol` is automatically
#'   scaled by `n / 1000` since the gradient is proportional to the
#'   sample size; for the random-intercept outer loop it is used as an
#'   absolute tolerance instead (see **Details**).
#' @param abstol maximum absolute change in parameter estimates from
#'   one iteration to the next before convergence can be declared; by
#'   default has no effect for the fixed-effects fit. For the
#'   random-intercept outer loop, this is capped at `1e-4` internally
#'   regardless of the value supplied, since the default of `1e10` was
#'   found to make this criterion a permanent no-op for that loop's
#'   only-linear convergence.
#' @param minstepsize used to specify when to abandon step-halving
#' @param tol singularity criterion. Default is typically `2e-16`
#' @param trace set to `TRUE` to print `-2` log likelihood,
#'   step-halving fraction, change in `-2` log likelihood, maximum
#'   absolute value of first derivative, and max absolute change in
#'   parameter estimates at each iteration. For the random-intercept
#'   outer loop, set to a value greater than `1` to also print
#'   per-outer-iteration mode-finding and AGQ diagnostics.
#' @param penalty.matrix a self-contained ready-to-use penalty matrix
#'   -- see `lrm`
#' @param weights a vector (same length as `y`) of possibly fractional
#'   case weights
#' @param normwt set to `TRUE` to scale `weights` so they sum to `n`,
#'   the length of `y`; useful for sample surveys as opposed to the
#'   default of frequency weighting
#' @param scale set to `TRUE` to subtract column means and divide by
#'   column standard deviations of `x` before fitting, and to
#'   back-solve for the un-normalized covariance matrix and regression
#'   coefficients. This can sometimes make the model converge for very
#'   large sample sizes where for example spline or polynomial
#'   component variables create scaling problems leading to loss of
#'   precision when accumulating sums of squares and crossproducts.
#' @param mscore set to `TRUE` to compute the sparse score matrix and
#'   store its elements as a list `mscore`. Not currently computed
#'   when `cluster` is given.
#' @param inclpen set to `FALSE` to not include the penalty matrix in
#'   the Hessian when the Hessian is being computed on transformed
#'   `x`, vs. adding the penalty after back-transforming. This should
#'   not matter.
#' @param y.precision when `y` is numeric, values may need to be
#'   rounded to avoid unpredictable behavior with `unique()` with
#'   floating-point numbers. Default is to 7 decimal places.
#' @param compstats set to `FALSE` to prevent the calculation of the
#'   vector of model statistics
#' @param onlydata set to `TRUE` to return the data used in model
#'   fitting as a list, without fitting the model
#' @param sigma.init starting value for the random-intercept standard
#'   deviation `sigma`, used only when `cluster` is given
#' @param maxit.outer maximum number of outer-loop iterations for the
#'   random-intercept fit (mode-finding / AGQ-reweighted Newton step /
#'   `sigma` profile search), used only when `cluster` is given
#' @param maxit.mode maximum number of Newton-Raphson iterations for
#'   the per-cluster random-intercept mode-finding inner loop, used
#'   only when `cluster` is given
#' @param nAGQ starting number of adaptive Gauss-Hermite quadrature
#'   nodes for the random-intercept fit, used only when `cluster` is
#'   given; escalated automatically as needed (see **Details**)
#' @param nAGQ.grid increasing sequence of candidate quadrature-node
#'   counts `nAGQ` is escalated through, used only when `cluster` is
#'   given
#' @param nAGQ.tol relative log-likelihood change, comparing the
#'   current `nAGQ` to the next candidate in `nAGQ.grid`, below which
#'   `nAGQ` escalation stops; used only when `cluster` is given
#' @param ... ignored
#'
#' @returns
#' a list with the following components, not counting all the
#' components produced by `orm.fit`:
#'
#' - `call`: calling expression
#' - `freq`: table of frequencies for `y` in order of increasing `y`
#' - `yunique`: vector of sorted unique values of `y`
#' - `stats`: vector with the following elements: number of
#'   observations used in the fit, effective sample size (`"ESS"`; see
#'   **Details** for a limitation when `cluster` is given), number of
#'   unique `y` values, median `y` from among the observations used in
#'   the fit, maximum absolute value of first derivative of log
#'   likelihood, model likelihood ratio chi-square, d.f., P-value,
#'   score chi-square and its P-value, Spearman's \eqn{\rho} rank
#'   correlation between linear predictor and `y` (if there is no
#'   censoring), Somers' \eqn{Dxy} rank correlation (if there is no
#'   censoring or only right censoring), the Nagelkerke \eqn{R^2}
#'   index, other \eqn{R^2} measures, the \eqn{g}-index, `gr` (the
#'   \eqn{g}-index on the ratio scale), and `pdm` (the mean absolute
#'   difference between 0.5 and the estimated probability that
#'   \eqn{y \geq} the marginal median). When `penalty.matrix` is
#'   present, the chi-square, d.f., and P-value are not corrected for
#'   the effective d.f. When `cluster` is given, the score-test-related
#'   elements are not currently available (see **Details**) and come
#'   back as `NA`.
#' - `fail`: set to `TRUE` if convergence failed (and `maxit>1`, or,
#'   when `cluster` is given, `maxit.outer>1`)
#' - `coefficients`: estimated parameters (intercepts then slopes)
#' - `family`, `famfunctions`: see [orm()]
#' - `deviance`: named vector of `-2` log likelihoods, in the order
#'   computed. Element names, in order (present or absent according to
#'   which of censoring/an offset/covariates/clustering actually
#'   apply -- see **Details** for why, when `cluster` is given and
#'   `p>0`, two additional, clustered elements are always computed):
#'   `"intercepts"` (always present); `"intercepts+offset"` (only if
#'   an offset variable is given); `"intercepts+x"` (only if `p>0`
#'   covariates are present); `"intercepts+random effects"` (only if
#'   `cluster` is given -- when `p>0` this is an extra, covariate-free
#'   clustered fit computed solely to give the next element a correct
#'   baseline; when `p==0` it is instead the final element, there
#'   being no `x` to add); `"intercepts+x+random effects"` (the final
#'   element whenever both `cluster` is given and `p>0`). The vector
#'   is correspondingly shorter whenever `cluster` is not given,
#'   ending at `"intercepts+x"` (or just `"intercepts"` if `p==0`).
#' - `lpe`: vector of per-observation likelihood probability elements.
#'   An observation's contribution to the log likelihood is the log of
#'   `lpe`. Not currently computed when `cluster` is given (returned
#'   as `NA`; see **Details**).
#' - `non.slopes`: number of intercepts in model
#' - `interceptRef`: the index of the middle (median) intercept used
#'   in computing the linear predictor and `var`
#' - `linear.predictors`: the linear predictor using the first
#'   intercept. When `cluster` is given, this is the fixed-effects-only
#'   linear predictor -- random effects are not added in, so that
#'   downstream fit statistics retain their ordinary interpretation.
#' - `penalty.matrix`: see above
#' - `info.matrix`: see [orm()] and [infoMxop()]. When `cluster` is
#'   given, this is the corrected (Louis's-identity-adjusted) sparse
#'   information matrix described in **Details**, with additional
#'   row(s)/column(s) in its `b`/`ab` pieces: a single `log(sigma)`
#'   element when `mre` is not given, or `log(sigma)` and `sigma2`
#'   together when it is. Use `infoMxop(fit$info.matrix,
#'   i='sigma_parameters')` to retrieve the variance of `log(sigma)`
#'   alone, or the `2x2` covariance matrix of `log(sigma)` and
#'   `sigma2` when both are present.
#' - `ncluster`: number of distinct clusters. Present only when
#'   `cluster` is given.
#' - `nAGQ`: the number of adaptive Gauss-Hermite quadrature nodes
#'   used for the final, converged fit (after any automatic
#'   escalation). Present only when `cluster` is given.
#' - `sigma`: the estimated random-intercept standard deviation.
#'   Present only when `cluster` is given *and* `mre` is not; when
#'   `mre` is given, `sigma1` and `sigma2` are reported instead (see
#'   below) and `sigma` is absent.
#' - `sigma1`, `sigma2`: the two estimated random-effect scale
#'   parameters described in **Details** (`sigma1` constrained
#'   positive, `sigma2` unconstrained). Present only when both
#'   `cluster` and `mre` are given.
#' - `gamma`: vector of estimated (posterior mode) random intercepts,
#'   one per cluster, in the order of the distinct levels of
#'   `cluster`. Present only when `cluster` is given.
#'
#' @author
#' Frank Harrell\cr
#' Department of Biostatistics, Vanderbilt University\cr
#' fh@fharrell.com
#'
#' @seealso
#' [orm()], [lrm()], [glm()], [gIndex()], [SparseM::solve()],
#' [recode2integer()], [infoMxop()]
#'
#' @examples
#' \dontrun{
#' # Fit an additive logistic model containing numeric predictors age,
#' # blood.pressure, and sex, assumed to be already properly coded and
#' # transformed
#' #
#' # fit <- orm.fit(cbind(age, blood.pressure, sex), death)
#'
#' # Fit a random-intercept version for repeated measurements on
#' # subject id
#' # fit <- orm.fit(cbind(age, blood.pressure, sex), death, cluster=id)
#'
#' # A Markov-1 model on repeated ordinal responses, where the linear
#' # predictor already includes a term for the previous response (e.g.
#' # via an added column of x, such as a spline in the lagged y). The
#' # recommended mre for this situation takes on only the values 0 and
#' # 1: 0 at each subject's FIRST follow-up time, 1 at every later one.
#' # This lets the fit determine, via the separately-estimated sigma2,
#' # how much the random effect should add to or subtract from the
#' # correlation the lag term alone induces, rather than assuming (as
#' # a plain, constantly-weighted random effect implicitly would) that
#' # its influence should be injected again at full, undiminished
#' # strength at every subsequent visit:
#' # mre <- ave(time, id, FUN = function(t) as.numeric(t > min(t)))
#' # fit <- orm.fit(cbind(prev.y, age, blood.pressure, sex), y,
#' #                 cluster=id, mre=mre)
#' }
#' @md
#' @export
#' @keywords models regression
#' @concept logistic regression model
orm.fit <- function(x=NULL, y, cluster=NULL,
                    family=c("logistic","probit","loglog","cloglog","cauchit"),
                    offset, initial,
                    opt_method=c('NR', 'LM'),
                    maxit=30L, eps=5e-4, gradtol=1e-3, abstol=1e10,
                    minstepsize=1e-2, tol=.Machine$double.eps, trace=FALSE,
                    penalty.matrix=NULL, weights=NULL, normwt=FALSE,
                    scale=FALSE, mscore=FALSE, inclpen=TRUE, y.precision = 7,
                    compstats=TRUE, onlydata=FALSE,
                    sigma.init=1.0, maxit.outer=100L, maxit.mode=30L,
                    nAGQ=7L, nAGQ.grid=c(7L, 11L, 15L, 21L, 31L, 45L, 63L),
                    nAGQ.tol=1e-5, mre=NULL, ...)
{
  cal        <- match.call()
  family     <- match.arg(family)
  opt_method <- match.arg(opt_method)

  debug2     <- getOption('orm.fit.debug2', FALSE)

  n <- NROW(y)
  if(! length(x)) {
      p     <- 0
      xname <- NULL
      x     <- 0.
    }
  else  {
      if(! is.matrix(x)) x <- as.matrix(x)
      dx <- dim(x)
      p  <- dx[2L]
      if(dx[1] != n) stop("x and y must have same number of rows")
      xname <- dimnames(x)[[2]]
      if(! length(xname)) xname <- paste("x[", 1 : p, "]", sep="")
  }
  ## Claude Sonnet 5 2026-09-08
  ## Diagnostic only: confirm p/x are already correct (or already
  ## broken) at the very top of orm.fit, before any of its own
  ## downstream logic runs.
  Hmisc::Fdebug('orm.re.debug')(p)
  Hmisc::Fdebug('orm.re.debug')(dim(x))

  len.penmat <- length(penalty.matrix)
  penpres    <- len.penmat && any(penalty.matrix != 0.)
  if(p == 0 && penpres) stop('may not specify penalty.matrix without predictors')
  if(penpres && any(dim(penalty.matrix) != p))
    stop(paste("penalty.matrix does not have", p, "rows and columns"))
  penmat <- if(! penpres) matrix(0e0, nrow=p, ncol=p) else penalty.matrix

  ## Extreme value type I dist = Gumbel maximum = exp(-exp(-x)) = MASS:::pgumbel
  ## Gumbel minimum = 1 - exp(-exp(x))
  families <- probabilityFamilies
  familiesDefined <- names(families)
  link  <- match(family, familiesDefined, nomatch=0)
  if(link == 0)
    stop('family must be one of ', paste(familiesDefined, collapse=' '))
  fam <- families[[family]]

  wtpres <- TRUE
  if(! length(weights)) {
    wtpres  <- FALSE
    normwt  <- FALSE
    weights <- rep(1.0, n)
  }
  if(length(weights) != n) stop('length of weights must equal length of y')
  if(normwt) weights <- weights * n / sum(weights)

  initial.there <- ! missing(initial) && length(initial)

  if(p > 0 && scale) {
    x      <- scale(x)
    scinfo <- attributes(x)[c('scaled:center', 'scaled:scale')]
    xbar   <- as.matrix(scinfo[[1]])
    xsd    <- as.matrix(scinfo[[2]])
    # Transform penalty matrix to new scale
    trans <- rbind(cbind(1., matrix(0., 1, p)),   # 1.: only dealing with middle intercept
                   cbind(- matrix(xbar / xsd, ncol=1),
                          diag(1. / as.vector(xsd), ncol=p)))
    if(penpres) penmat <- t(trans[-1, -1]) %*% penmat %*% trans[-1, -1]
    }


  # model.extract which calls model.response may not keep Ocens class
  isOcens <- NCOL(y) == 2
  yupper  <- NULL
  npsurv  <- NULL
  ctype   <- integer(n)

  if(isOcens) {
    Y       <- y                # save original values before Ocens2ord
    aY      <- attributes(Y)
    y       <- Ocens2ord(Y)
    a       <- attributes(y)
    ylabel  <- aY$label
    uni     <- aY$units
    ylevels <- a$levels
    yupper  <- a$upper
    npsurv  <- a$npsurv
    if(! initial.there && ! length(npsurv))
      stop('npsurv must be present when not specifying initial to orm.fit')
    Ncens1  <- a$Ncens1
    Ncens2  <- a$Ncens2
    # Claude Sonnet 5 High 2026-07-10          2 lines
    rt_cens_beyond <- a$rt_cens_beyond   # deprecated; see tailinfo
    tailinfo       <- a$tail
    k       <- length(ylevels) - 1
    YO      <- extractCodedOcens(y, what=4, ivalues=TRUE)   # ivalues=TRUE -> a, b are [0, k]
    ctype   <- YO$ctype
    y       <- YO$a
    y2      <- YO$b
    # Ocens uses -Inf and Inf for left and right censoring; later we transform to integer codes -1, k+1
    y2[ctype == 2] <- k + 1
    y [ctype == 1] <- -1
    numy    <- a$freq
    mediany <- a$median
    kmid    <- which.min(abs(ylevels[-1] - mediany))
    anycens <- any(ctype > 0)
    ranges  <- a$ranges
  } else {  # regular variable
    ylabel    <- label(y)
    uni       <- units(y)
    mul       <- 10^y.precision
    yrange    <- if(is.numeric(y)) range(round(y * mul)) / mul   # be consistent with recode2integer
    ranges    <- if(is.numeric(y)) list(y=yrange, u=yrange, c=c(NA, NA))
    w         <- recode2integer(y, precision=y.precision)
    y         <- w$y - 1
    y2        <- y
    ylevels   <- w$ylevels
    k         <- length(ylevels) - 1
    kmid      <- max(w$whichmedian - 1L, 1L)
    numy      <- w$freq
    mediany   <- w$median
    Ncens1    <- Ncens2 <- c(left=0L, right=0L, interval=0L)
    anycens   <- FALSE
    rt_cens_beyond <- NULL
    # Claude Sonnet 5 High 2026-07-10          1 line
    tailinfo       <- NULL
  }

  intcens <- 1L * any(ctype == 3)

  if(k == 1) kmid <- 1

  ylev <- ylevels[-1L]
  if(length(yupper)) ylev <- ifelse(yupper[-1] > ylev, paste(ylev, '-', yupper[-1]), ylev)
  iname <- if(k == 1) "Intercept" else paste0("y>=", ylev)
  # Claude Sonnet 5 High 2026-07-10          4 lines
  # An added level for right censoring at or beyond the highest uncensored
  # value represents P(Y > highest uncensored value); label it accordingly
  if(k > 1 && length(tailinfo) && tailinfo$type == 'right')
    iname[tailinfo$index - 1L] <- paste0('y>', ylevels[tailinfo$index - 1L])
  name  <- c(iname, xname)

  if(onlydata) return(
    if(isOcens) list(Y=cbind(y, y2), Ncens1=Ncens1, Ncens2=Ncens2, k=k, iname=iname)
    else list(y=y, k=k, iname=iname) )

  if(missing(offset) || ! length(offset) || (length(offset) == 1 && offset == 0.))
    offset <- rep(0., n)
  ofpres <- ! all(offset == 0.)
  if(ofpres && length(offset) != n) stop("offset and y must have same length")

  if(n < 3) stop("must have >=3 non-missing observations")

  nv <- p + k

  # These are used for initial intercept estimates when there is no censoring (OK)
  # and also in R2 statistics
  sumwty <- tapply(weights, y, sum)
  sumwt  <- sum(sumwty)
  if(anycens) {
    ess <- n - sum(Ncens1)
    Nw  <- n
  } else {
    rfrq   <- sumwty / sumwt
    ess    <- n * (1 - sum(rfrq ^ 3))
    Nw     <- sumwt
  }

  finverse <- eval(fam[2])

  if(! initial.there) {
    if(anycens) {
      # If only censored values are right censored, then MLEs of intercepts are exactly
      # the link function of Kaplan-Meier estimates.  We're ignoring weights.
      # This should also work for only left censoring, and reasonable results
      # are expected for interval and mixed censoring using a Turnbull-type estimator
      # Claude Sonnet 5 High 2026-07-10          2 lines
      # computed by Ocens2ord's internal self-consistency EM (Ocens_npmle) and
      # stored in the npsurv object; the legacy cons='data' path uses icenReg::ic_np instead.
      pp <- npsurv$surv[-1]
    } else {
      ncum   <- rev(cumsum(rev(sumwty)))[2 : (k + 1)]
      pp     <- ncum / sumwt
    }
    initial <- finverse(pp)
    if(ofpres) initial <- initial - mean(offset)
    initial <- c(initial, rep(0., p))
 } else if(is.list(initial) && length(initial) == 2) {
   # initial is an npsurv object from a previous run, which is used when bootstrapping
   # Use linear interpolation/extrapolation from the stored times to the current times
   interp_initial <- approxExtrap(initial$time, initial$surv, xout=ylevels)$y
   if(debug2) {
     du <- cbind(time    = initial$time, surv       = initial$surv,
                 ylevels = ylevels,      interpsurv = interp_initial)
     saveRDS(du, '/tmp/du.rds')
     prn(du[1:min(20, length(ylevels)),], 'orm.fit')
   }
   initial <- finverse(interp_initial[-1])
   if(ofpres) initial <- initial - mean(offset)
   initial <- c(initial, rep(0., p))
 }

  ## Claude Sonnet 5 2026-09-07
  ## Every element pushed into loglik (-> the returned `deviance`
  ## vector) is named at the point it is computed, regardless of which
  ## combination of censoring/offset/covariates/clustering branches
  ## actually run -- see orm.fit's own @returns documentation for the
  ## full, ordered list of possible names.
  if(! anycens) { loglik <- -2 * sum(sumwty * log(sumwty / sum(sumwty)))
                  names(loglik) <- 'intercepts' }

  if(anycens || (p==0 & ! ofpres)) {
    z <- ormfit(NULL, y, y2, k, intcens, initial=initial[1 : k],
                offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method, maxit=maxit,
                tolsolve=tol, objtol=eps, gradtol=gradtol, paramtol=abstol,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    kof    <- z$coef
    if(debug2) {
      prn(intcens)
      ww <- cbind(initial=initial[1 : k], kof)
      print(utils::head(ww, 10)); print(utils::tail(ww, 10))
    }
    loglik <- z$loglik
    names(loglik) <- 'intercepts'
    initial <- c(kof, rep(0., p))
    info   <- z$info
  }

  if(ofpres) {
    ## Fit model with only intercept(s) and offset
    ## Check that lrm.fit uses penmat in this context   ??
    z <- ormfit(NULL, y, y2, k, intcens, initial=initial[1 : k],
                offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method, maxit=maxit,
                tolsolve=tol, objtol=eps, gradtol=gradtol, paramtol=abstol,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    kof    <- z$coef
    loglik <- c(loglik, z$loglik)
    names(loglik)[length(loglik)] <- 'intercepts+offset'
    initial <- c(z$coef, rep(0., p))
    if(p == 0) info <- z$info
  }

  if(p > 0) {
    # Fit model with intercept(s), offset, covariables
    z <- ormfit(x, y, y2, k, intcens, initial=initial, offset=offset, wt=weights,
                penmat=penmat, opt_method=opt_method,
                maxit=maxit, tolsolve=tol, objtol=eps,
                gradtol=gradtol, paramtol=abstol, mscore=mscore,
                trace=trace, link=link, iname=iname, xname=xname)
    if(z$fail) return(structure(list(fail=TRUE), class="orm"))
    loglik <- c(loglik, z$loglik)
    names(loglik)[length(loglik)] <- 'intercepts+x'
    kof  <- z$coef
    info <- z$info
    # Compute linear predictor before unscaling beta, as x is scaled
    lp <- matxv(x, kof, kint=kmid)

    if(scale) {
      betas         <- kof[- (1 : k)]
      kof[1 : k]    <- kof[1 : k] - sum(betas * xbar / xsd)
      kof[-(1 : k)] <- betas / xsd
      xbar          <- as.vector(xbar)
      names(xbar)   <- xname
      xsd           <- as.vector(xsd)
      names(xsd)    <- xname
      info$scale    <- list(mean=xbar, sd=xsd)
    }
  } else lp <- rep(kof[kmid], n)

  ## Claude Sonnet 5 2026-09-08
  ## Diagnostic only: pinpoint whether p/x/kof are already broken
  ## HERE, at the boundary between orm.fit's own fixed-effects logic
  ## above and the cluster dispatch below, or whether they were still
  ## correct up to this point and something in the dispatch block
  ## itself is at fault.
  deb <- Hmisc::Fdebug('orm.re.debug')
  deb(p); deb(dim(x)); deb(length(kof)); deb(kof)

  ## Claude Sonnet 5 2026-08-30          entire block
  ## Random-intercept dispatch. If cluster is present, override the
  ## fixed-effects-only z/kof/info/lp computed above with a call to
  ## ormrfit, using the ALREADY-refined fixed-effects fit (kof, from
  ## whichever combination of the censoring/offset/covariate branches
  ## above actually ran) as ormrfit's warm start -- this reuses ALL of
  ## orm.fit's existing initial-value-refinement machinery
  ## unconditionally (its cost is paid once whether or not clusters
  ## end up being used) and eliminates the double front-end-data-
  ## preparation pass a separate orm.rfit() calling a nested orm.fit()
  ## previously required. ormfit/ormrfit remain the only places actual
  ## model-fitting math happens; this is purely dispatch.
  if(length(cluster)) {
    cluster.factor <- as.factor(cluster)
    cluster.int    <- as.integer(cluster.factor)
    cluster.levels <- levels(cluster.factor)
    nc             <- length(cluster.levels)
    if(length(cluster.int) != n) stop('cluster must have the same length as y')

    ## Claude Sonnet 5 2026-08-30          4 lines
    ## Do not pass abstol straight through as ormrfit's paramtol:
    ## abstol's default (1e10) is correct/preserved for ormfit (its
    ## quadratic Newton convergence means this loose default never
    ## actually binds), but ormrfit's block-alternating convergence
    ## needs paramtol to actually be tight to mean anything (see the
    ## earlier "convergence too loose" fix) -- cap it here so a user
    ## who never touched abstol still gets that fix, while a user who
    ## explicitly passed a tighter abstol has it respected for both.
    ## Claude Sonnet 5 2026-09-07          entire block
    ## For a correct beta-specific LR test under clustering, also fit
    ## an INTERCEPT-ONLY (p=0) clustered model -- the SAME random-
    ## effects structure (mre, if present), just without x -- to get a
    ## clustered null deviance actually comparable to the clustered
    ## full model's. Appending BOTH the clustered null and full
    ## deviances (rather than just the full one, the previous
    ## behavior) keeps the stats block's existing "the last two
    ## elements of loglik are [null, full]" logic correct without
    ## needing to touch it: previously, with only the full clustered
    ## deviance appended, loglik[length-1] silently grabbed the
    ## NON-clustered full-model deviance instead of any null deviance
    ## at all, so the reported "Model L.R." was actually testing "does
    ## the random effect improve fit over the fixed-effects-only
    ## model" rather than "is beta significant" -- confirmed directly:
    ## it bore no relationship to beta's own Wald chi-square.
    if(p > 0) {
      ## Claude Sonnet 5 2026-09-08          entire block
      ## Refit a proper, converged non-clustered INTERCEPT-ONLY model
      ## first, rather than warm-starting the covariate-free clustered
      ## null model directly from kof[1:k] -- the full model's own
      ## intercepts, jointly optimized WITH p covariates. Removing
      ## those covariates while keeping the same intercepts is often a
      ## poor starting point: with no beta contribution left to
      ## explain any of the spread across categories, some adjacent
      ## intercepts can end up implying a near-zero (or, under
      ## floating-point rounding, exactly zero) category probability
      ## that was perfectly fine under the full model's own linear
      ## predictor. Confirmed directly to matter in practice under
      ## probit specifically: its faster-saturating CDF turns this
      ## mismatch into repeated, though individually recoverable
      ## (step-halved) salloc=999 hits in the clustered null-model fit.
      ##
      ## Claude Sonnet 5 2026-09-08          9 lines
      ## kof[1:k] is NOT used as this refit's OWN starting point either
      ## -- confirmed directly (from a real failure) that it can
      ## already be invalid there too, with zero NR iterations even
      ## attempted, meaning ordinary step-halving cannot help: it only
      ## rescues an overshooting STEP, not a starting point that is
      ## itself already outside the valid region. finverse(pp) -- the
      ## same link-appropriate, purely marginal-frequency-based
      ## intercept-only starting guess orm.fit itself computes at the
      ## very top for its OWN first fit, independent of any covariate
      ## fit -- is recomputed fresh here (rather than relying on
      ## whatever `initial` happened to be if the user supplied their
      ## own) as a properly-scaled starting point instead.
      pp0 <- if(anycens) npsurv$surv[-1]
             else { ncum0 <- rev(cumsum(rev(sumwty)))[2 : (k + 1)]; ncum0 / sumwt }
      init_null <- finverse(pp0)
      if(ofpres) init_null <- init_null - mean(offset)
      z_null_fixed <- ormfit(NULL, y, y2, k, intcens, initial=init_null,
                             offset=offset, wt=weights, penmat=matrix(0., 0, 0),
                             opt_method=opt_method, maxit=maxit,
                             tolsolve=tol, objtol=eps, gradtol=gradtol,
                             paramtol=abstol, trace=trace, link=link,
                             iname=iname, xname=character(0))
      if(z_null_fixed$fail) return(structure(list(fail=TRUE), class="orm"))
      zr0 <- ormrfit(x=matrix(0., n, 0), y=y, y2=y2, k=k, intcens=intcens,
                    cluster=cluster.int, nc=nc,
                    initial=z_null_fixed$coef, sigma.init=sigma.init,
                    offset=offset, wt=weights, penmat=matrix(0., 0, 0),
                    maxit=maxit, maxit.outer=maxit.outer, maxit.mode=maxit.mode,
                    objtol=eps, gradtol=gradtol, paramtol=min(abstol, 1e-4),
                    tolsolve=tol, minstepsize=minstepsize, trace=trace,
                    link=link, iname=iname, xname=character(0),
                    nAGQ=nAGQ, nAGQ.grid=nAGQ.grid, nAGQ.tol=nAGQ.tol, mre=mre)
      if(zr0$fail) return(structure(list(fail=TRUE), class="orm"))
      loglik <- c(loglik, zr0$loglik)
      names(loglik)[length(loglik)] <- 'intercepts+random effects'
    }

    zr <- ormrfit(x=x, y=y, y2=y2, k=k, intcens=intcens,
                 cluster=cluster.int, nc=nc,
                 initial=kof, sigma.init=sigma.init,
                 offset=offset, wt=weights, penmat=penmat,
                 maxit=maxit, maxit.outer=maxit.outer, maxit.mode=maxit.mode,
                 objtol=eps, gradtol=gradtol, paramtol=min(abstol, 1e-4),
                 tolsolve=tol, minstepsize=minstepsize, trace=trace,
                 link=link, iname=iname, xname=xname,
                 nAGQ=nAGQ, nAGQ.grid=nAGQ.grid, nAGQ.tol=nAGQ.tol, mre=mre)
    if(zr$fail) return(structure(list(fail=TRUE), class="orm"))
    z      <- zr
    kof    <- zr$coef
    info   <- zr$info
    loglik <- c(loglik, zr$loglik)
    names(loglik)[length(loglik)] <- if(p > 0) 'intercepts+x+random effects' else 'intercepts+random effects'
    ## Stats block below treats random effects as absent: lp is the
    ## fixed-effects-only linear predictor, computed via the exact
    ## same formula used above for the non-clustered case, not
    ## anything gamma-adjusted.
    lp     <- if(p > 0) matxv(x, kof, kint=kmid) else rep(kof[kmid], n)
  } else if(length(mre))
    stop('mre requires cluster to also be specified')

  # Add second derivative of penalty function if needed, on the original scale
  if(! inclpen && penpres)
    info$b <- info$b - penalty.matrix

  names(kof) <- name

  stats <- NULL
  if(compstats) {
    if(p == 0) {llnull <- loglik[length(loglik)]; model.lr <- 0e0}
    else {
      llnull   <- loglik[length(loglik) - 1L]
      model.lr <- llnull - loglik[length(loglik)]
    }
    model.df <- p
    if(initial.there || maxit == 1)
      model.p <- score <- score.p <- NA
    else {
      score <- z$score
      if(model.df == 0)
        model.p <- score.p <- 1.
      else {
        model.p <- 1. - pchisq(model.lr, model.df)
        score.p <- 1. - pchisq(score,    model.df)
      }
    }

    # Effective sample size ess did not count censored observations
    # They do contain some information, especially for interval-censored ones
    # with small intervals.  Adjust ess for partial contributions.

    ## Claude Sonnet 5 2026-08-30          5 lines
    ## z$lpe is currently always NA for a clustered fit (ormrfit does
    ## not yet compute a per-observation marginal probability -- a
    ## known, flagged gap, not an oversight). Previously this let NA
    ## propagate into ess and then into R2Measures(), which errors
    ## outright on an NA input rather than returning NA gracefully --
    ## crashing EVERY clustered+censored fit, not just leaving one
    ## stat missing as intended. Skip the adjustment (leaving ess at
    ## its already-computed, pre-adjustment value) when lpe isn't
    ## actually available, rather than corrupting ess and crashing.
    if(anycens && length(cluster))
      warning('ESS adjustment for censored observations is not yet ',
              'available for clustered (random-intercept) fits -- ess ',
              'reflects only the uncensored-observation count.')
    if(anycens && ! length(cluster)) {
      nuncens <- ess   # based on original counts, not Ocens modifications
      ll      <- -2e0 * log(z$lpe)
      a       <- sum(ll[ctype == 0])   # ctype references censoring after Ocens adjustments ?? TODO
      # Compute multiplier that makes ll for uncensored obs sum to the number uncensored
      b       <- nuncens / sum(ll[ctype == 0])
      # Compute this scaled contribution of censored obs
      essc    <- b * sum(ll[ctype > 0])
      ess     <- ess + essc
    }

    r2     <- 1. - exp(- model.lr / Nw)
    r2.max <- 1. - exp(- llnull / Nw)
    r2     <- r2 / r2.max
    r2m    <- R2Measures(model.lr, model.df, Nw, ess)
    if(k > 1L) attr(lp, 'intercepts') <- kmid
    g  <- GiniMd(lp)
    ## compute average |difference| between 0.5 and the condition
    ## probability of being >= marginal median
    cump <- eval(fam[1])
    pdm <- mean(abs(cump(lp) - 0.5))
    rho <- if(anycens) NA else if(p == 0 || diff(range(lp)) == 0e0) 0e0 else cor(rank(lp), rank(y))
    ## Somewhat faster:
    ## rho <- .Fortran('rcorr', cbind(lp, y), as.integer(n), 2L, 2L, r=double(4),
    ##                 integer(4), double(n), double(n), double(n), double(n),
    ##                 double(n), integer(n), PACKAGE='Hmisc')$r[2]

    cindex <- NA
    if(! anycens || (Ncens1['left'] == 0 && Ncens1['interval'] == 0))
      cindex <- suppressWarnings(concordancefit(if(anycens) Ocens2Surv(Y) else y, lp)$concordance)
    Dxy    <- 2 * (cindex - .5)

    stats <- c(n, ess, length(numy), mediany, z$dmax, model.lr, model.df,
              model.p, score, score.p, rho, Dxy, r2, r2m, g, exp(g), pdm)

    nam <- c("Obs", "ESS", "Distinct Y", "Median Y", "Max Deriv",
            "Model L.R.", "d.f.", "P", "Score", "Score P",
            "rho", "Dxy", "R2", names(r2m), "g", "gr", "pdm")
    names(stats) <- nam
    }

  info$iname <- iname
  if(! length(cluster)) info$xname <- xname

  retlist <- list(call              = cal,
                  freq              = numy,
                  yunique           = ylevels,
                  yupper            = yupper,             # NULL if no censoring that produced yupper > ylevels
                  ylabel            = ylabel,
                  units             = uni,
                  Ncens1            = if(isOcens) Ncens1,
                  Ncens2            = if(isOcens) Ncens2,
                  # n.risk          = if(any(ctype > 0)) n.risk,
                  ranges            = ranges,             # ranges attribute from Ocens or computed above
                  # Claude Sonnet 5 High 2026-07-10          2 lines
                  rt_cens_beyond    = rt_cens_beyond,     # deprecated; see tail
                  tail              = tailinfo,
                  stats             = stats,
                  coefficients      = kof,
                  var               = NULL,
                  u                 = z$u,
                  lpe               = z$lpe,
                  iter              = z$iter,
                  family            = family,
                  famfunctions      = fam,
                  deviance          = loglik,
                  non.slopes        = k,
                  interceptRef      = kmid,
                  linear.predictors = structure(lp, intercepts=kmid),
                  penalty.matrix    = if(penpres) penalty.matrix,
                  weights           = if(wtpres) weights,
                  xbar              = if(p > 0 && scale) xbar,
                  xsd               = if(p > 0 && scale) xsd,
                  info.matrix       = info,
                  mscore            = z$mscore,
                  fail              = FALSE)

  ## Claude Sonnet 5 2026-08-30          4 lines
  ## Added only when cluster is present -- no placeholders otherwise,
  ## matching every other cluster-only field (info.matrix's own shape
  ## already varies this way for the interval-censoring case).
  ## Claude Sonnet 5 2026-09-06          2 lines
  ## zr$sigma is NULL when mre was used (sigma1/sigma2 replace it
  ## entirely -- see ormrfit's own return value) -- attach whichever
  ## of the two is actually present, never a placeholder for the other.
  if(length(cluster)) {
    retlist$ncluster <- nc
    retlist$nAGQ     <- zr$nAGQ
    retlist$gamma    <- zr$gamma
    if(length(mre)) {
      retlist$sigma1 <- zr$sigma1
      retlist$sigma2 <- zr$sigma2
    } else
      retlist$sigma  <- zr$sigma
  }

  class(retlist) <- 'orm'
  retlist
}

ormfit <-
  function(x, y, y2, k, intcens, link, initial,
           offset=rep(0., n), wt=rep(1., n), penmat=matrix(0., p, p), opt_method='NR',
           maxit=30L, objtol=5e-4, gradtol=1e-3, paramtol=1e10, tolsolve=.Machine$double.eps,
           minstepsize=1e-2, mscore=FALSE, trace=FALSE, iname, xname) {

# y  =  -1 for left censored observation
# y2 = k+1 for right censored observation
# Uncensored observations have y = y2 = 0, 1, ..., k
# intcens = 1 if there are any interval censored observations

if(getOption('orm.fit.debug', FALSE)) try <- function(x) x

n <- length(y)
p <- length(initial) - k

if(k > 1 && any(diff(initial[1:k]) >= 0))
  stop('initial values for intercepts are not in descending order')

storage.mode(x)       <- 'double'
storage.mode(y)       <- 'integer'
storage.mode(y2)      <- 'integer'
storage.mode(k)       <- 'integer'
storage.mode(intcens) <- 'integer'
storage.mode(p)       <- 'integer'
storage.mode(initial) <- 'double'
storage.mode(offset)  <- 'double'
storage.mode(wt)      <- 'double'
storage.mode(penmat)  <- 'double'
storage.mode(link)    <- 'integer'

# Claude Sonnet 5 2026-08-30          9 lines
# Precompute intercept-index bookkeeping once per fit (not once per
# Newton-Raphson iteration) -- pure function of y, y2, k. See ormidx.f90.
idx <- .Fortran(F_ormidx, n, k, y, y2,
                 ia=integer(n), ia2=integer(n), sgn=numeric(n),
                 ib=integer(n), nb=integer(1), salloc=integer(1))
if(idx$salloc != 0)
  stop('Censoring values encountered that are not handled (ormidx code ',
       idx$salloc, ')')
nb <- idx$nb
ib <- idx$ib[1 : nb]

if(getOption('rms.fit.debug', FALSE)) try <- function(x) x

rfort <- function(theta, what=3L, mscore=FALSE,
                  debug=as.integer(getOption('orm.fit.debug', 0L)),
                  debug2=getOption('orm.fit.debug2', FALSE)) {
  p <- as.integer(length(theta) - k)
  nai <- as.integer(if(intcens) 1000000 else 0)
  if(debug2) prn(c(k, length(theta), p), 'rfort')
  if(debug) {
    a <- llist(n, k, p, y, y2, link, intcens, nai, what, debug)
    s <- sapply(a, storage.mode)
    if(any(s != 'integer')) stop(s)
    ac <- llist(x, offset, wt, penmat, theta[1:k], theta[-(1:k)],
                logL=numeric(1))
    sc <- sapply(ac, storage.mode)
    if(any(sc != 'double')) stop(s)
    g <- function(x) if(is.matrix(x)) paste(dim(x), collapse='x') else length(x)
    print(sapply(c(a, ac), g), quote=FALSE)
  }
  nu <- if(mscore) n * (2L + p) else 0L
  # Claude Sonnet 5 2026-08-30          2 lines
  # ia, ia2, sgn, ib, nb now passed through (computed once, above) rather
  # than being recomputed inside ormll on every call. y, y2 are no longer
  # passed at all -- ormll's only uses of them were debug-diagnostic
  # (now served by ia/ia2/sgn instead) and a dead assignment.
  w <- .Fortran(F_ormll, n, k, p, x, idx$ia, idx$ia2, idx$sgn, ib, nb,
                offset, wt, penmat,
                link=link, theta[1:k], theta[-(1:k)],
                logL=numeric(1), grad=numeric(k + p), lpe=numeric(n),
                a=matrix(0e0, (1 - intcens) * k, 2), b=matrix(0e0, p, p), ab=matrix(0e0, k, p),
                intcens, row=integer(nai), col=integer(nai), ai=numeric(nai), nai=nai, ne=integer(1),
                urow=integer(nu), ucol=integer(nu), um=numeric(nu), nu=nu, nuu=integer(1),
                what=what, debug=as.integer(debug), 1L, salloc=integer(1))
  if(debug) prn(w$salloc)

  if(w$salloc == 999) {   # zero or negative probability in likelihood calculation
    w$logL = Inf          # triggers step-halving in MLE updating loop
    return(w)
  }
  if(w$salloc == 998)
    stop('Censoring values encountered that are not handled')
  if(w$salloc == 997)
    stop('More than 1,000,000 elements needed in the intercepts part of the information matrix due\n',
         'to the variety of interval censored observations')
  if(w$salloc == 996)
    stop('more than ', nu, ' elements needed in the score matrix')
  if(w$salloc != 0)
    stop('Failed dynamic array allocation in Fortran subroutine ormll: code ', w$salloc)
  if(intcens && what == 3L) {
    ne    <- w$ne
    w$a   <- list(row = w$row[1 : ne], col = w$col[1 : ne], a = w$ai[1 : ne])
    w$row <- w$col <- w$ai <- w$nai <- w$ne <- NULL
    }
  if(mscore) {
    nuu      <- w$nuu
    w$mscore <- Matrix::sparseMatrix(w$urow[1 : nuu], w$ucol[1 : nuu], x = w$um[1 : nuu], dims=c(n, k + p))
  }
  w
  }

  if(missing(x) || ! length(x) || p == 0) {
    x <- 0.
    p <- 0L
  }

  nv <- k + p
  if(length(initial) < nv)
    initial <- c(initial, rep(0., nv - length(initial)))
  if(trace > 2) prn(initial)

m <- function(x) max(abs(x))

if(maxit == 1) {
  w      <- rfort(initial)
  # Information matrix is negative Hessian on LL scale
  if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
  info <- list(a = w$a, b = - w$b, ab = - w$ab,
               iname=iname, xname=xname)

  res <- list(coefficients = initial,
              loglik       = w$logL,
              info         = info,
              u            = w$grad,
              lpe          = w$lpe,
              dmax=m(w$grad), score=NA,
              iter=1, fail=FALSE, class='orm')
    return(res)
  }

  theta      <- initial # Initialize the parameter vector
  oldobj     <- 1e12
  score.test <- NA

  gradtol <- gradtol * n / 1e3

  # Newton-Raphson MLE with step-halving, initial draft generated by ChatGPT
  if(opt_method == 'NR') {
    for (iter in 1:maxit) {
      w <- rfort(theta)
      if(iter == 1) objf <- w$logL
      gradient <- w$grad
      hess     <- infoMxop(w[c('a', 'b', 'ab')])

      # Newton-Raphson step

      delta <- try(Matrix::solve(hess, gradient, tol=tolsolve))
      # Runs amazingly slow if Matrix:: is omitted; prob. not using Matrix
      if(inherits(delta, 'try-error')) {
        message('singular Hessian matrix')
        return(list(fail=TRUE))
      }

      if(trace > 0)
        cat('NR iteration:', iter, '  -2LL:', format(objf, nsmall=4),
            '  Max |gradient|:', m(gradient),
            '  Max |change in parameters|:', m(delta), '\n', sep='')

      if(opt_method == 'NR' && is.na(score.test) && p > 0 &&
        all(theta[- (1 : k)] == 0.))
        score.test <- - gradient %*% delta

      step_size <- 1.0           # Initialize step size for step-halving

      # Step-halving loop
      while (TRUE) {
        new_theta <- theta - step_size * delta # Update parameter vector
        wd <- which(diff(new_theta[1 : k]) >= 0e0)
        if(length(wd)) {
          if(trace > 0) cat('new_theta out of order for intercepts',
              paste(wd + 1, collapse=' '), 'forced step-halving\n')
          objfnew <- Inf
        }
        else objfnew   <- rfort(new_theta, what=1L)$logL
        if(trace > 1)
          cat('Old, new, old - new -2 LL:', objf, objfnew, objf - objfnew, '\n')
        if (! is.finite(objfnew) || objfnew > objf + objtol / 10.) {
          # Objective function failed to be reduced or is infinite
          step_size <- step_size / 2e0         # Reduce the step size
          if(trace > 0) cat('Step size reduced to', step_size, '\n')
          if(step_size < minstepsize) {
            message('Step size ', step_size, ' has reduced below minstepsize=',
                    minstepsize,
                    ' without improving log likelihood; fitting stopped')
            return(list(fail=TRUE))
          }
        } else {
          theta  <- new_theta                   # Accept the new parameter vector
          oldobj <- objf
          objf   <- objfnew
          if(trace > 2) prn(theta)
          break
        }
      }

      # Convergence check - must meet 3 criteria
      if((objf <= oldobj + objtol / 10. && (oldobj - objf < objtol)) &&
        (m(gradient) < gradtol) &&
        (m(delta)    < paramtol)) {
        # Compute final information matrix (in 3 parts) since not computed
        # since Newton-Raphson updating
        w <- rfort(theta, mscore=mscore)
        if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
        info <- list(a = w$a, b = - w$b, ab = - w$ab,
                    iname=iname, xname=xname)

        return(list(coef           = theta,
                    loglik         = w$logL,
                    u              = w$grad,
                    mscore         = w$mscore,
                    lpe            = w$lpe,
                    info           = info,
                    objchange      = oldobj - w$logL,
                    dmax           = m(w$grad),
                    maxparamchange = m(delta),
                    score          = score.test,
                    iter           = iter,
                    fail           = FALSE) )
        }
    }

    msg <- paste('Reached', maxit, 'iterations without convergence\nChange in -2LL:',
      oldobj -objf, ' Max |gradient|:', m(gradient),
      ' Max |change in parameters|:', m(delta))
    message(msg)
    return(list(fail=TRUE))

  } else {    # L-M

  lambda   <- 1e-3    # hard-wired for L-M
  oldobj   <- 1e12
  objf     <- NA      # needed in case no H_damped is ever positive definite
  w        <- rfort(theta)
  gradient <- w$grad
  H        <- infoMxop(w[c('a', 'b', 'ab')])

  for (iter in 1:maxit) {
    H_damped <- H + lambda * Matrix::Diagonal(x = Matrix::diag(H))
    delta    <- try(Matrix::solve(H_damped, gradient, tol=tolsolve))
    if(inherits(delta, 'try-error')) {
      # Increase lambda if Hessian is ill-conditioned
      lambda <- lambda * 10.
      next
      }

    theta_new <- theta - delta
    objf      <- rfort(theta_new, what=1L)$logL
    if(trace > 0)
      cat('LM iteration:', iter, '  -2LL:', format(objf, nsmall=4),
          '  Max |gradient|:', m(gradient),
          '  Max |change in parameters|:', m(delta), '\n', sep='')
    if(trace > 1)
      cat('Old, new, old - new -2 LL:', oldobj, objf, oldobj - objf, '\n')

    if(is.finite(objf) &&
       (objf <= oldobj + objtol / 10. && (oldobj - objf < objtol)) &&
       (m(gradient) < gradtol) &&
       (m(delta)    < paramtol)) break

    if(is.finite(objf) && (objf < oldobj)) {
      # Accept the step and decrease lambda
      # Claude Sonnet 5 2026-08-31
      # Was lambda / 10 (symmetric with the rejection-side * 10 below).
      # Confirmed by direct reproduction on a real, moderately ill-
      # conditioned (but NOT saddle-point) dataset that the symmetric
      # rule can lock into a permanent limit cycle: dividing by 10
      # immediately after a success jumps straight from "just enough
      # damping to succeed" to "too little, overshoots" with no
      # intermediate ground, so lambda bounces between the same two
      # values forever (confirmed directly: iterations 12-40 cycling
      # lambda between 1e-4 and 1e-5, -2LL frozen, never escaping).
      # Asymmetric adjustment (slow to relax damping, fast to
      # reintroduce it on a rejection) breaks this: verified on the
      # same dataset to converge cleanly in 19 iterations, -2LL
      # decreasing smoothly throughout, max|delta| shrinking
      # monotonically to 0 -- no remaining cycle.
      theta    <- theta_new
      oldobj   <- objf
      w        <- rfort(theta)
      gradient <- w$grad
      H        <- infoMxop(w[c('a', 'b', 'ab')])
      lambda   <- lambda / 1.5
    } else {
      # Reject the step and increase lambda
      lambda <- lambda * 10.
    }
  }
  if(iter == maxit) {
    msg <- paste('Reached', maxit, 'iterations without convergence\n-2LL:',
      objf, ' Max |gradient|:', m(gradient))
    message(msg)
    return(list(fail=TRUE))
  }
  w     <- rfort(theta, mscore=mscore)
  if(intcens) w$a$a <- - w$a$a else w$a <- - w$a
  info  <- list(a = w$a, b = - w$b, ab = - w$ab,
                iname=iname, xname=xname)

  return(list(coef           = theta,
              loglik         = w$logL,
              u              = w$grad,
              mscore         = w$mscore,
              lpe            = w$lpe,
              info           = info,
              objchange      = objf - w$logL,
              dmax           = m(w$grad),
              maxparamchange = m(delta),
              score          = NA,
              iter           = iter,
              fail           = FALSE) )
  }   # End M-L
}

## Note: deriv and deriv2 below are no longer used as are hard-coded into ormll
## Expressions are used because if using a function that calls plogis(),
## the .C call for plogis will can result in R losing the environment
## of the C code.
## The 5 expression elements are cumprob, inverse, deriv, deriv2, and
## deriv as a function only of x

## Extreme value type I dist = Gumbel maximum = exp(-exp(-x)) = MASS:::pgumbel
## Gumbel minimum = 1 - exp(-exp(x))
probabilityFamilies <-
  list(logistic =
         expression(function(x) plogis(x),
                    function(x) qlogis(x),
                    function(x, f) f*(1-f),
                    function(x, f, deriv) f*(1-3*f+2*f*f),
                    function(x) {f <- plogis(x); f*(1-f)}),
       probit =
         expression(function(x) pnorm(x),
                    function(x) qnorm(x),
                    function(x, f) dnorm(x),
                    function(x, f, deriv) - deriv * x,
                    function(x) dnorm(x)),
       loglog =
         expression(function(x) exp(-exp(-x)),
                    function(x) -log(-log(x)),
                    function(x, f) exp(-x-exp(-x)),
                    function(x, f, deriv)
                      ifelse(abs(x) > 200, 0,
                             exp(-x - exp(-x)) * (-1 + exp(-x))),
                    function(x) exp(-x-exp(-x))),
       cloglog =
         expression(function(x) 1-exp(-exp(x)),
                    function(x) log(-log(1-x)),
                    function(x, f) exp(x-exp(x)),
                    function(x, f, deriv)
                      ifelse(abs(x) > 200, 0, deriv * ( 1 - exp( x))) ,
                    function(x) exp(x-exp(x))),
       cauchit =
         expression(function(x) pcauchy(x),
                    function(x) qcauchy(x),
                    function(x, f) dcauchy(x),
                    function(x, f, deriv) -2 * x * ((1 + x*x)^(-2)) / pi,
                    function(x) dcauchy(x))
  )

## Check:
## P(x) = plogis(x); P'(x) = P(x) - P(x)^2
## d <- function(x) plogis(x) - 3*plogis(x)^2 + 2*plogis(x)^3
## x <- seq(-3, 3, length.out=150)
## plot(x, d(x), type='l')
## ad <- c(NA,diff(dlogis(x))/(x[2]-x[1]))
## lines(x, ad, col='red')
