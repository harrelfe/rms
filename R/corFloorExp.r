#' Floor-plus-exponential-decay correlation structure for `nlme`
#'
#' @description
#' An [nlme::corStruct] class for within-subject correlations that decay
#' exponentially from 1 at zero lag down to a fixed, non-negative floor as
#' lag grows, with the decay rate itself allowed to shift with how far
#' along in follow-up a pair of observations sits (non-isotropy).
#'
#' @details
#' For two observation times \eqn{t_1, t_2} on the same subject, with
#' \eqn{L = |t_1-t_2|} and \eqn{m = (t_1+t_2)/2}, the modeled correlation is
#'
#' \deqn{
#'   \mathrm{corr}(t_1, t_2) \;=\;
#'   k \;+\; (1-k)\,\exp\!\bigl(-\mathrm{rate}(m)\, L\bigr),
#'   \qquad \mathrm{rate}(m) = \exp(a + b\,m)
#' }{
#'   corr(t1, t2) = k + (1-k) * exp( -rate(m) * L ),
#'   rate(m) = exp(a + b*m)
#' }
#'
#' This is the classical "shared random effect plus serial decay" mixture
#' (a continuous-time analogue of a random intercept superimposed on an
#' Ornstein-Uhlenbeck / exponential-decay process): `k` is the correlation
#' contributed by a persistent, subject-level component that never fades,
#' and `(1-k)` is the share of correlation carried by a component that
#' decays with lag at rate `rate(m)`.
#'
#' **Boundary conditions hold by construction, not just at the fit.**
#' At \eqn{L=0}, \eqn{\exp(-\mathrm{rate}(m)\cdot 0)=1}{exp(-rate(m)*0) = 1}, so
#' \eqn{\mathrm{corr}=k+(1-k)=1}{corr = k + (1-k) = 1} for every `m`  -- 
#' same-time observations are always perfectly correlated. As
#' \eqn{L\to\infty}{L -> Inf},
#' \eqn{\exp(-\mathrm{rate}(m)\,L)\to 0}{exp(-rate(m)*L) -> 0}, so
#' \eqn{\mathrm{corr}\to k}{corr -> k} for every `m`  --  `k` is a literal
#' floor, not an artifact of the fitted range.
#'
#' **Parametrization and validity.** All three parameters are estimated on
#' an unconstrained (real-line) scale, matching what `nlme`'s optimizer
#' wants internally:
#' \itemize{
#'   \item `k` is stored on the logit scale; the reported floor is
#'     \eqn{k = \mathrm{plogis}(k^{*})}{k = plogis(k*)}, which is
#'     automatically in \eqn{(0,1)}  --  non-negative, and strictly below 1
#'     as required for a valid floor.
#'   \item `a` sets the baseline decay rate at \eqn{m=0} via
#'     \eqn{\mathrm{rate}=\exp(a)}{rate = exp(a)}; the exponential link
#'     keeps the rate positive for any `a`, so the decaying component can
#'     never invert direction or overshoot.
#'   \item `b` is the non-isotropy term: it lets \eqn{\mathrm{rate}(m)}{rate(m)}
#'     grow or shrink with the pair's mean time `m` rather than staying
#'     fixed across the study.
#' }
#' Every entry of the resulting matrix is confined to \eqn{[k, 1]} by
#' construction, so  --  unlike a raw linear-in-covariates form  --  this
#' structure cannot produce an out-of-range correlation for any parameter
#' values.
#'
#' **Positive-definiteness.** The isotropic sub-case (`b=0`) is provably
#' valid for *any* set of time points: it is a convex combination of two
#' known-valid correlation matrices (the all-ones matrix, from the
#' persistent shared component, and the standard exponential/OU kernel
#' from the decaying component), and a convex combination of valid
#' correlation matrices is always itself a valid correlation matrix. Once
#' `b != 0`, the decay rate differs pair by pair and that guarantee no
#' longer strictly applies, so  --  as a numerical safeguard rather than a
#' modeling choice  --  the per-subject matrix is checked at every
#' evaluation and, if not positive definite, is replaced by the nearest
#' positive definite, unit-diagonal correlation matrix (negative
#' eigenvalues clipped to a small positive floor, reconstructed, and
#' rescaled to a unit diagonal). In simulation this guard is essentially
#' never invoked at realistic parameter magnitudes; it exists to keep the
#' optimizer from crashing if a trial step wanders into an aggressively
#' non-isotropic region.
#'
#' As with any `corStruct`, all Cholesky-based factorization and
#' log-determinant machinery needed by `gls()`/`lme()` is supplied by
#' `nlme`'s generic `corFactor.corStruct` default operating on the matrix
#' this class returns.
#'
#' **Combining with a random intercept.** If this structure is used
#' inside [nlme::lme()] together with `random = ~ 1 | subject`, the
#' explicit random intercept and the floor `k` both represent a
#' persistent, lag-independent source of within-subject correlation, so
#' they compete for the same signal; expect `k` to shrink relative to a
#' `gls()`-only fit, and treat the two as jointly, not separately,
#' identified in that case.
#'
#' To test whether the non-isotropy term is needed (\eqn{b=0}), fit the
#' nested 2-parameter version with `b` fixed at 0 and compare by
#' likelihood ratio.
#'
#' @param value Numeric vector of length 3, `c(k, a, b)`, giving the
#'   starting values on the scales described above (`k` on the logit
#'   scale, `a` and `b` already unconstrained). Defaults to
#'   `c(k = 0, a = 0, b = 0)`, i.e. a starting floor of 0.5 and a starting
#'   baseline decay rate of 1.
#' @param form A one-sided formula of the form `~ time | subject`
#'   identifying the time covariate and, optionally, the grouping
#'   (subject) factor, following the usual [nlme::corStruct] convention.
#' @param fixed Logical. If `TRUE`, `k`, `a`, and `b` are held at `value`
#'   rather than estimated.
#'
#' @return An (uninitialized) correlation structure object of class
#'   `c("corFloorExp", "corStruct")`, suitable for the `correlation`
#'   argument of [nlme::gls()] or [nlme::lme()].
#'
#' @seealso [nlme::corClasses], [nlme::gls()], [nlme::lme()],
#'   [nlme::corExp()], [nlme::corCompSymm()]
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' times <- replicate(60, sort(sample(0:12, sample(4:7, 1))), simplify = FALSE)
#' dat <- do.call(rbind, Map(function(tt, id) {
#'   data.frame(subject = id, time = tt,
#'              y = 2 + 0.1 * tt + rnorm(length(tt)))
#' }, times, seq_along(times)))
#' dat$subject <- factor(dat$subject)
#'
#' fit <- nlme::gls(y ~ time, data = dat,
#'                   correlation = corFloorExp(form = ~ time | subject))
#' coef(fit$modelStruct$corStruct, unconstrained = FALSE)
#' }
#'
#' @export
corFloorExp <- function(value = c(k = 0, a = 0, b = 0), form = ~ time | subject,
                         fixed = FALSE) {
  if (length(value) != 3)
    stop("'value' must have 3 elements: k (floor, logit scale), a (log baseline rate), b (non-isotropy)")
  value <- as.numeric(value)
  names(value) <- c("k", "a", "b")
  attr(value, "formula") <- form
  attr(value, "fixed")   <- fixed
  class(value) <- c("corFloorExp", "corStruct")
  value
}

## ------------------------------------------------------------------
## Internal machinery below: not exported. corMatrix.corFloorExp supplies
## the correlation matrix; nlme's generic corFactor.corStruct default
## handles the Cholesky factorization/log-determinant from there.
## ------------------------------------------------------------------

## Project a symmetric matrix onto the nearest valid PD, unit-diagonal
## correlation matrix; used as a guard when non-isotropy (b != 0) pushes
## the raw formula outside the provably-valid isotropic case (see Details).
#' @noRd
.nearestCorPD <- function(R, eps = 1e-6) {
  ee <- eigen((R + t(R)) / 2, symmetric = TRUE)
  if (min(ee$values) > eps) return(R)
  vals <- pmax(ee$values, eps)
  Rn <- ee$vectors %*% (vals * t(ee$vectors))
  d  <- sqrt(diag(Rn))
  Rn <- Rn / outer(d, d)
  diag(Rn) <- 1
  (Rn + t(Rn)) / 2
}

#' @exportS3Method nlme::Initialize
#' @noRd
Initialize.corFloorExp <- function(object, data, ...) {
  object <- NextMethod()
  object
}

#' @exportS3Method nlme::corMatrix
#' @noRd
corMatrix.corFloorExp <- function(object, covariate = getCovariate(object),
                                   corr = TRUE, ...) {
  par <- as.vector(object)
  kstar <- par[1]; a <- par[2]; b <- par[3]
  k <- plogis(kstar)

  buildOne <- function(tt) {
    n <- length(tt)
    if (n <= 1) return(matrix(1, n, n))
    lag  <- abs(outer(tt, tt, "-"))
    mid  <- outer(tt, tt, "+") / 2
    rate <- exp(a + b * mid)
    R <- k + (1 - k) * exp(-rate * lag)
    diag(R) <- 1
    .nearestCorPD(R)
  }

  if (is.list(covariate)) {
    val <- lapply(covariate, buildOne)
    if (length(val) == 1) val <- val[[1]] else names(val) <- names(covariate)
  } else {
    val <- buildOne(covariate)
  }
  val
}

#' @exportS3Method base::coef
#' @noRd
coef.corFloorExp <- function(object, unconstrained = TRUE, ...) {
  if (attr(object, "fixed") && unconstrained) return(numeric(0))
  val <- as.vector(object)
  if (!unconstrained) val[1] <- plogis(val[1])
  names(val) <- c("k", "a", "b")
  val
}

#' @exportS3Method base::coef<-
#' @noRd
"coef<-.corFloorExp" <- function(object, ..., value) {
  if (length(value) != 3)
    stop("cannot change the length of the parameter of a \"corFloorExp\" object")
  object[] <- value
  object
}

#' @exportS3Method base::print
#' @noRd
print.corFloorExp <- function(x, ...) {
  if (length(aux <- coef(x, unconstrained = FALSE)) > 0) {
    cat("Correlation structure of class corFloorExp representing\n")
    cat("  corr(t1,t2) = k + (1-k)*exp(-exp(a+b*(t1+t2)/2)*|t1-t2|)\n")
    print(aux, ...)
  } else {
    cat("Uninitialized correlation structure of class corFloorExp\n")
  }
  invisible(x)
}
