##' Precompute Intercept Subscripts for Ordinal Regression Likelihoods
##'
##' Calls the compiled routine `ormidx` (the same one `orm.fit` uses) to
##' compute, once per fit, the per-observation intercept subscripts that enter
##' the log-likelihood of a cumulative probability ordinal regression model.
##' The result depends only on `y`, `y2`, and `k`, not on any parameter value,
##' so it can be computed once and then passed to other code, for example as
##' data to a Stan program in the `rmsb` package, so that the likelihood does
##' not need conditionals on the type of censoring.
##'
##' `y` and `y2` use the integer coding produced by [Ocens2ord()] and used
##' internally by [orm.fit()], where `k` is the number of intercepts, i.e.,
##' one less than the number of distinct ordinal levels:
##' * `y == y2` in `0:k`: uncensored observation
##' * `y == -1`: left censored, with upper bound `y2`
##' * `y2 == k + 1`: right censored, with lower bound `y`
##' * `0 <= y < y2 <= k`: interval censored
##'
##' The contribution of observation `i` to the likelihood is
##' `p1 - p2` where `p1 = F(alpha[ia[i]] + Xb)` if `sgn[i] == 1` or
##' `1 - F(alpha[ia[i]] + Xb)` if `sgn[i] == -1`, and `p2 = 1 - F(alpha[ia2[i]] + Xb)`
##' for the observations in `ib` (the only ones with `ia2 > 0`), where `F` is the
##' cumulative distribution function of the link family.  All subscripts are
##' 1-based so they can be used directly in R or in Stan.
##' @param y integer-valued vector (see Details), no `NA`s
##' @param y2 integer-valued vector of the same length as `y`, defaulting to `y` (no censoring)
##' @param k number of intercepts (a single positive integer)
##' @return a list with integer vectors `ia` and `ia2` (length `n`; `ia2` is 0 when there is
##' no second intercept), numeric vector `sgn` (length `n`, `1` or `-1`), integer vector `ib`
##' (observation numbers for which `ia2 > 0`, of length `nb`), and integer `nb`.
##' @seealso [orm.fit()], [Ocens2ord()]
##' @md
##' @author Frank Harrell
##' @examples
##' # 5 ordinal levels => k=4 intercepts; last 2 observations are right and left censored
##' y  <- c(0L, 2L, 4L, 3L, -1L)
##' y2 <- c(0L, 2L, 4L, 5L,  2L)
##' ormIdx(y, y2, k=4L)
##' @export
ormIdx <- function(y, y2 = y, k) {
  n <- length(y)
  if(length(y2) != n) stop('y and y2 must have the same length')
  if(anyNA(y) || anyNA(y2)) stop('y and y2 may not contain NA')
  if(length(k) != 1L || is.na(k) || k < 1 || k != round(k))
    stop('k must be a single positive integer')
  yi  <- as.integer(y)
  y2i <- as.integer(y2)
  if(any(yi != y) || any(y2i != y2)) stop('y and y2 must be integer-valued')
  idx <- .Fortran(F_ormidx, as.integer(n), as.integer(k), yi, y2i,
                  ia=integer(n), ia2=integer(n), sgn=numeric(n),
                  ib=integer(n), nb=integer(1), salloc=integer(1))
  if(idx$salloc != 0)
    stop('Censoring values encountered that are not handled (ormidx code ',
         idx$salloc, ')')
  list(ia=idx$ia, ia2=idx$ia2, sgn=idx$sgn,
       ib=idx$ib[seq_len(idx$nb)], nb=idx$nb)
}
