##' Censored Ordinal Variable
##'
##' Creates a 2-column integer matrix that handles left- right- and interval-censored ordinal or continuous values for use in [rmsb::blrm()] and [orm()].  A pair of values `[a, b]` represents an interval-censored value known to be in the interval `[a, b]` inclusive of `a` and `b`.  It is assumed that all distinct values are observed as uncensored for at least one observation.  When both input variables are `factor`s it is assumed that the one with the higher number of levels is the one that correctly specifies the order of levels, and that the other variable does not contain any additional levels.  If the variables are not `factor`s it is assumed their original values provide the orderings.  A left-censored point is is coded as having `-Inf` as a lower limit, and a right-censored point is coded as having `Inf` as an upper limit.   As with most censored-data methods, modeling functions assumes that censoring is independent of the response variable values that would have been measured had censoring not occurred.  `Ocens` creates a 2-column integer matrix suitable for ordinal regression.  Attributes of the returned object give more information.
##'
##' The upper limit of a left-censored observation is coded as the minimum uncensored value that is at or above the original upper limit.
##' The lower limit of a right-censored observation is coded as the maximum uncensored value that is at or below the original lower limit.
##' For interval-censored values, the limits are from uncensored values at the limits or just interior to them.  These steps are taken
##' because there are no intercept parameters in ordinal models for censored points that do not also correspond to uncensored values.
##'
##' @param a vector representing a `factor`, numeric, or alphabetically ordered character strings.  Censoring points have values of `-Inf`.
##' @param b like `a`.  If omitted, it copies `a`, representing nothing but uncensored values.  Censoring points have values of `Inf`.
##' @param precision when `a` and `b` are numeric, values may need to be rounded to avoid unpredictable behavior with \code{unique()} with floating-point numbers. Default is to 7 decimal places.
##' @verbose set to `FALSE` to suppress information messages
##' @return a 2-column integer matrix of class `"Ocens"` with an attribute `levels` (ordered).  When the original variables were `factor`s, these are factor levels, otherwise are numerically or alphabetically sorted distinct (over `a` and `b` combined) values.  When the variables are not factors and are numeric, another attributes `median` and `range` are also returned.  `median` is the median of the uncensored values on the origiinal scale.  `range` is a 2-vector range of original data values before adjustments to closest uncensored points.  When the variables are factor or character, the median of the integer versions of variables for uncensored observations is returned as attribute `mid`.  A final attribute `freq` is the vector of frequencies of occurrences of all uncensored values.  `freq` aligns with `levels`.
##' @author Frank Harrell
##' @export
Ocens <- function(a, b=a, precision=7, verbose=TRUE) {
  nf <- is.factor(a) + is.factor(b)
  if(nf == 1)
    stop('one of a and b is a factor and the other is not')

  if(nf == 0 && is.numeric(a) && is.numeric(b)) {
    mul <- 1e0
    z <- c(a, b);  z <- z[is.finite(z)]
    if(any(z %% 1 != 0)) {   # see recode2integer
      a <- round(a * 10^precision)
      b <- round(b * 10^precision)
      mul <- 10^-precision
    }
    uncensored <- ! is.na(a) & ! is.na(b) & (a == b)
    if(! any(uncensored)) stop('no uncensored observations')

    if(any(is.infinite(a) & is.infinite(b))) stop('an observation has infinite values for both values')
  
    # Set censoring points to next uncensored values
    u <- sort(unique(a[! is.na(a + b) & a == b]))    # distinct uncensored values
    A <- a; B <- b
    # Left censored: [-Inf, b] -> [-Inf, next highest u]
    for(j in which(is.infinite(a))) {
      jj <- which(u >= b[j])
      if(! length(jj)) stop('Left censoring point ', mul*b[j],
                            ' has no uncensored observations above it')
      y <- min(u[jj])
      if(b[j] != y) b[j] <- y
    }
    # Right censored: [a, max] -> [next lowest u, max]
    for(j in which(is.infinite(b))) {
      jj <- which(u <= a[j])
      if(! length(jj)) stop('Right censoring point ', mul*a[j],
                            ' has no uncensored observations below it')
      y <- max(u[jj])
      if(a[j] != y) a[j] <- y
    }
    # Interval censored: [a, b] -> [next lowest u, next highest u]
    for(j in which(is.finite(a) & is.finite(b) & a != b)) {
      if(! any(u <= a[j]) || ! any(u >= b[j]))
        stop('no uncensored observations outside intervals of interval censored values',
             ' [', mul*a[j], ',', mul*b[j], ']')
      yl <- max(u[u <= a[j]])
      yu <- min(u[u >= b[j]])
      if(a[j] != yl || b[j] != yu) {
        a[j] <- yl
        b[j] <- yu
      }
    }

    w <- 1 * (A != a) + 2 * (B != b)
    if(verbose && any(w > 0)) {
      d <- data.frame(Obs=1 : length(a),
                      'Old Lower'=A*mul, 'Old Upper'=B*mul,
                      'New Lower'=a*mul, 'New Upper'=b*mul,
                      check.names=FALSE)
      if(any(w == 1)) {
        cat('Modified lower values to next lower uncensored value:\n\n')
        print(subset(d, w == 1, c('Obs', 'Old Lower', 'New Lower')), row.names=FALSE)
      }
      if(any(w == 2)) {
        cat('Modified upper values to next upper uncensored value:\n\n')
        print(subset(d, w == 2, c('Obs', 'Old Upper', 'New Upper')), row.names=FALSE)
      }
      if(any(w == 3)) {
        cat('Modified lower and upper values to closest uncensored values outside the interval:\n\n')
        print(subset(d, w == 3), row.names=FALSE)
      }
    }

    # Since neither variable is a factor we can assume they are ordered
    # numerics.  Create an nx2 matrix of integers
    if(any(b < a)) stop('some values of b are less than corresponding a values')
    ymed <- median(a[uncensored], na.rm=TRUE) * mul
    ai    <- match(a, u, nomatch=0)
    bi    <- match(b, u, nomatch=0)
    if(any((ai == 0) != (a == -Inf))) stop('program logic error 1 in Ocens')
    if(any((bi == 0) != (b ==  Inf))) stop('program logic error 2 in Ocens')
    ai[ai == 0] <- -Inf
    bi[bi == 0] <-  Inf
    freq <- tabulate(ai[uncensored], nbins=length(u))
    y <- cbind(ai, bi)
    dimnames(y) <- list(NULL, NULL)
    return(structure(y,
                     class  = 'Ocens',
                     levels = u * mul,
                     freq   = freq,
                     median = ymed,
                     range  = range(z)))
  }

  uncensored <- ! is.na(a) & ! is.na(b) & (a == b)
  if(! any(uncensored)) stop('no uncensored observations')

  alev <- levels(a)
  blev <- levels(b)
  ## Cannot just pool the levels because ordering would not be preserved
  if(length(alev) >= length(blev)) {
    master <- alev
    other  <- blev
  } else{
    master <- blev
    other  <- alev
  }
  if(any(other %nin% master))
    stop('a variable has a level not found in the other variable')
  a <- match(as.character(a), master)
  b <- match(as.character(b), master)
  if(any(b < a)) stop('some values of b are less than corresponding a values')
  freq <- tabulate(a[uncensored], nbins=length(master))
  mid  <- quantile(a[uncensored], probs=.5, type=1L)
  structure(cbind(a, b), class='Ocens', levels=master, freq=freq, mid=mid)
}

##' Convert `Ocens` Object to Data Frame to Facilitate Subset
##'
##' Converts an `Ocens` object to a data frame so that subsetting will preserve all needed attributes
##' @param x an `Ocens` object
##' @param row.names optional vector of row names
##' @param optional set to `TRUE` if needed
##' @param ... ignored
##' @return data frame containing a 2-column integer matrix with attributes
##' @author Frank Harrell
##' @export
as.data.frame.Ocens <- function(x, row.names = NULL, optional = FALSE, ...) {
  nrows <- NROW(x)
  row.names <- if(optional) character(nrows) else as.character(1:nrows)
  value <- list(x)
  if(! optional) names(value) <- deparse(substitute(x))[[1]]
  structure(value, row.names=row.names, class='data.frame')
}

##' Subset Method for `Ocens` Objects
##'
##' Subsets an `Ocens` object, preserving its special attributes.  Attributes are not updated.  In the future such updating should be implemented.
##' @title Ocens
##' @param x an `Ocens` object
##' @param rows logical or integer vector
##' @param cols logical or integer vector
##' @param ... ignored
##' @return new `Ocens` object
##' @author Frank Harrell
##' @md
##' @export
'[.Ocens' <- function(x, rows=1:d[1], cols=1:d[2], ...) {
  d <- dim(x)
  at <- attributes(x)[c('levels', 'freq', 'median', 'range', 'mid')]
  x <- NextMethod('[')
  attributes(x) <- c(attributes(x), at)
  class(x) <- 'Ocens'
  x
  }

Ocens2Surv <- function(Y) {
  y  <- Y[, 1]
  y2 <- Y[, 2]
  su <- survival::Surv
  if(all(y == y2)) return(su(y))  # no censoring
  i <- which(is.finite(y) & is.finite(y2))
  w <- 1 * any(is.infinite(y)) + 2 * any(is.infinite(y2)) + 4 * any(y[i] != y2[i])
  # Only left censoring:
  if(w == 1)      su(y2, event=y == y2, type='left')
  else if(w == 2) su(y,  event=y == y2, type='right')
  else if(w == 4) su(y, event=rep(3, length(y)), time2=y2,       type='interval')
  else            su(y, time2=y2,       type='interval2')
}

# g <- function(a, b) {
#   s <- Ocens2Surv(cbind(a, b))
#   print(s)
#   km.quick(s, interval='>=')
# }
# g(1:3, 1:3)
# g(c(-Inf, 2, 3), c(2.5, 2, 3))
# g(1:3, c(1, 2, Inf))
# g(c(1, 4, 7), c(2, 4, 8))
# g(c(-Inf, 2,4, 6), c(3, 3, 4, Inf))
# a <- c(-Inf, 2, 1, 4, 3)
# b <- c(   3, 3, 1, 5, 3)
# g(a, b)
