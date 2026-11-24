#' Title survest.orm
#'
#' @param fit result of `orm`
#' @param newdata data frame defining covariate settings
#' @param linear.predictors linear predictor vector using the reference intercept
#' @param x design matrix
#' @param times times for which estimates are desired; defaults to estimating probabilities of T > t for all uncensored times
#' @param fun optional transformation of survival probabilities
#' @param loglog set to `TRUE` to use the log-log transformatino
#' @param conf.int a number between 0-1 with the default of 0.95; set to 0 to not compute CLs
#' @param what specify `what='parallel'` to compute the survival probability at the observed linear predictor and time values, both varying; all possible combinations of these are then not created
#' @param ... ignored
#'
#' @returns a data frame with variables `time, surv`.  If `conf.int > 0` the data also contains `lower, upper`.  The variable `Xrow` indicates the row of the design matrix or the linear predictor element used in getting the current data frame row estimate.
#' @export
#'
#' @md
#' @author Frank Harrell
#' @examples
#' # See survest.psm
survest.orm <- function(fit, newdata=NULL, linear.predictors=NULL, x=NULL, times=NULL,
                        fun, loglog=FALSE, conf.int=0.95,
                        what=c("survival", "parallel"),
                        ...) {

  what <- match.arg(what)
  if(what=='parallel') conf.int <- FALSE

  S <- Survival(fit)

  # Claude Sonnet 5 2026-11-24          entire block
  # A yde-fitted orm() object's Survival() returns a closure with a
  # DIFFERENT signature (times, newdata, conf.int, ...) than the
  # ordinary one (times, lp, X, ...) -- see ormYdeSurvival() in
  # predict.lrm.s; a Y-dependent effect cannot be represented by a
  # single linear predictor value, so that closure needs the actual
  # covariate data frame. The code below (written only for the
  # ordinary closure) ends with S(times, linear.predictors, X=x, ...);
  # linear.predictors is never actually set from x (it stays at its
  # NULL default throughout this function), so for a yde fit that NULL
  # was being passed POSITIONALLY into the yde closure's own newdata
  # argument -- confirmed directly (reported by Frank from
  # survplot(g, trt); reproduced the exact error, then isolated it to
  # this positional mismatch) -- giving a confusing "newdata must be
  # given" error even when newdata WAS supplied to survest() itself
  # (it just never reached the closure). Route directly to the yde
  # closure's own newdata= argument instead, bypassing all of the
  # ordinary closure's lp/X-oriented logic below, none of which
  # applies to a yde fit.
  if(length(fit$qyde) && fit$qyde > 0) {
    if(what == 'parallel')
      stop('what="parallel" is not yet implemented for an orm() fit using yde/ytrans')
    if(! length(newdata))
      stop('survest() for an orm() fit using yde/ytrans must be called with ',
           'newdata (the original covariate data frame) -- a bare linear ',
           'predictor or design matrix (linear.predictors=/x=) cannot ',
           'represent a Y-dependent effect, which varies by level of Y.')
    return(S(newdata=newdata, times=times, conf.int=conf.int, forcedf=TRUE, zero=TRUE))
  }

  if(missing(fun)) fun <- if(loglog) function(x) logb(ifelse(x == 0 | x == 1, NA, x))
  else function(x) x

  if(conf.int > 0 && ! missing(linear.predictors)) {
    warning('conf.int set to 0 since linear.predictors specified')
    conf.int <- 0
  }

  p <- length(fit$coef) - num.intercepts(fit)

  lx <- length(x)
  ln <- length(newdata)
  lt <- length(times)

  if(what == 'parallel' && ! lt) stop('times must be given when what="parallel"')

  if(! length(linear.predictors)) {
    linear.predictors <- if(p == 0) 0
    if(p > 0 && ! lx && ! ln) {
      if(what != 'parallel' && ! lt)
        stop('specify times= if using linear predictors from fit')
      linear.predictors <- fit$linear.predictors
      if(conf.int > 0)
        stop("may not specify conf.int > 0 unless x or newdata given")
      rnam <- names(linear.predictors)
    }
    else {
      if(! lx) x <- predict(fit, newdata, type="x")
      rnam <- dimnames(x)[[1]]
    }
  }
  else  rnam <- names(linear.predictors)

  if(what == 'parallel') 
    return(S(times, linear.predictors, parallel=TRUE))

  n <- length(linear.predictors)

  if(n > 1 & ! lt)
    warning("should specify times if getting predictions for >1 obs.")

  S(times, linear.predictors, X=x, conf.int=conf.int, forcedf=TRUE, zero=TRUE)
  }
