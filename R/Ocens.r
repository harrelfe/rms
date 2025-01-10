##' Censored Ordinal Variable
##'
##' Creates a 2-column integer matrix that handles left- right- and interval-censored ordinal or continuous values for use in [rmsb::blrm()] and in the future [orm()].  A pair of values `[a, b]` represents an interval-censored value known to be in the interval `[a, b]` inclusive of `a` and `b`.  It is assumed that all distinct values are observed as uncensored for at least one observation.  When both input variables are `factor`s it is assume that the one with the higher number of levels is the one that correctly specifies the order of levels, and that the other variable does not contain any additional levels.  If the variables are not `factor`s it is assumed their original values provide the orderings.  Since all values that form the left or right endpoints of an interval censored value must be represented in the data, a left-censored point is is coded as `a=1` and a right-censored point is coded as `b` equal to the maximum observed value.  If the maximum observed value is not really the maximum possible value, everything still works except that predictions involving values above the highest observed value cannot be made.  As with most censored-data methods, modeling functions assumes that censoring is independent of the response variable values that would have been measured had censoring not occurred.
##' @param a vector representing a `factor`, numeric, or alphabetically ordered character strings
##' @param b like `a`.  If omitted, it copies `a`, representing nothing but uncensored values
##' @param precision when `a` and `b` are numeric, values may need to be rounded to avoid unpredictable behavior with \code{unique()} with floating-point numbers. Default is to 7 decimal places.
##' @return a 2-column integer matrix of class `"Ocens"` with an attribute `levels` (ordered).  When the original variables were `factor`s, these are factor levels, otherwise are numerically or alphabetically sorted distinct (over `a` and `b` combined) values.  When the variables are not factors and are numeric, another attribute `median` is also returned.  This is the median of the uncensored values.  When the variables are factor or character, the median of the integer versions of variables for uncensored observations is returned as attribute `mid`.  A final attribute `freq` is the vector of frequencies of occurrences of all uncensored values.  `freq` aligns with `levels`.
##' @author Frank Harrell
##' @export
Ocens <- function(a, b=a, precision=7) {
  nf <- is.factor(a) + is.factor(b)
  if(nf == 1)
    stop('one of a and b is a factor and the other is not')

  uncensored <- ! is.na(a) & ! is.na(b) & (a == b)
  if(! any(uncensored)) stop('no uncensored observations')
  ymed <- if(is.numeric(a)) median(a[uncensored])

  if(nf == 0) {
    num <- is.numeric(a) + is.numeric(b)
    if(num == 1) stop('one of a and b is numeric and the other is not')
    if((num == 2) && (any(c(a, b) %% 1 != 0))) {   # see recode2integer
      a <- round(a * 10^precision) * 10^-precision
      b <- round(b * 10^precision) * 10^-precision
    }
    ## Since neither variable is a factor we can assume they are ordered
    ## numerics or use standard character string ordering
    ## Combine distinct values and create an nx2 matrix of integers
    u <- sort(unique(c(a, b)))
    u <- u[! is.na(u)]
    if(any(b < a)) stop('some values of b are less than corresponding a values')
    a    <- match(a, u)
    b    <- match(b, u)
    freq <- tabulate(a[uncensored], nbins=length(u))
    return(structure(cbind(a, b), class='Ocens', levels=u, freq=freq, median=ymed))
  }
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
  at <- attributes(x)[c('levels', 'freq', 'mid')]
  x <- NextMethod('[')
  attributes(x) <- c(attributes(x), at)
  class(x) <- 'Ocens'
  x
  }
