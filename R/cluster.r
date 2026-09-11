##' Cluster Function for Random Effects
##'
##' Used by `orm` and `blrm` to signal a categorical variable to generate random effects.
##' @title cluster
##' @param x a vector representing a categorical variable
##' @return `x` unchanged
##' @author Frank Harrell
##' @md
##' @export
cluster <- function(x) x

##' Mixing Function for Dual Random Effect Sigmas
##'
##' Used by `orm` to specify weights for mixing two scaled normal(0,1) random effects.  The main purpose of `mre` is to get a more realistic correlation structure in the raw data space when fitting Markov ordinal models with lags in the raw Y values.  What works well in practice is to set `mre` to 0.0 for an observation measured at the earliest follow-up time, and `mre` to `1.0` when the observation is taken at any time after this initial visit time.
##' @title mix_re
##' @param mre mixing proportion for dual sigmas for random effects with [orm()].  These are typically a function of time in a Markov model, scaled to be between 0 and 1 inclusive.  When `mre` is zero the multiplier of normal(0,1) random effects is the first sigma, and when `mre` is 1.0 the multiplier is the second sigma.  Otherwise it is a weighted combination of sigmas.  The first sigma is constrained to be positive.
##' @return `mre` unchanged
##' @author Frank Harrell
##' @md
##' @export
##' @examples
##' # mix_re(ifelse(t == 1, 0, 1))
mix_re <- function(mre) mre
