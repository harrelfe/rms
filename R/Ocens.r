##' Censored Ordinal Variable
##'
##' Creates a 2-column integer matrix that handles left- right- and interval-censored ordinal or continuous values for use in [rmsb::blrm()] and [orm()].  A pair of values `[a, b]` represents an interval-censored value known to be in the interval `[a, b]` inclusive of `a` and `b`. Left censored values are coded as `(-Infinity, b)` and right-censored as `(a, Infinity)`, both of these intervals being open at the finite endpoints. Open left and right censoring intervals are created by adding a small increment (subtracting for left censoring) to `a` or `b`. When this occurs at the outer limits, new ordinal categories will be created by `orm` to capture the real and unique information in outer censored values.  For example if the highest uncensored value is 10 and there is a right-censored value in the data at 10, a new category `10+` is created, separate from the category for `10`.  So it is assumed that if an exact value of 10 was observed, the pair of values for that observation would not be coded as `(10, Infinity)`.
##'
##' The intervals that drive the coding of the input data into numeric ordinal levels are the Turnbull intervals computed by the non-exported `findMaximalIntersections` function in the `icenReg` package, which handles all three types of censoring.  These are defined in the `levels` and `upper` attributes of the object returned by `Ocens`.  Sometimes consecutive Turnbull intervals contain the same statistical information likelihood function-wise, leading to the same survival estimates over two ore more consecutive intervals.  This leads to zero probabilities of involved ordinal values, preventing `orm` from computing a valid log-likeliihood.  A limited about of interval consolidation is done by `Ocens` to alleviate this problem.  Depending on the value of `cons` this consolidation is done by intervals (preferred) or by changing the raw data.  If `verbose=TRUE`, information about the actions taken is printed.
##'
##' When both input variables are `factor`s it is assumed that the one with the higher number of levels is the one that correctly specifies the order of levels, and that the other variable does not contain any additional levels.  If the variables are not `factor`s it is assumed their original values provide the orderings.  A left-censored point is is coded as having `-Inf` as a lower limit, and a right-censored point is coded as having `Inf` as an upper limit.   As with most censored-data methods, modeling functions assumes that censoring is independent of the response variable values that would have been measured had censoring not occurred.  `Ocens` creates a 2-column integer matrix suitable for ordinal regression.  Attributes of the returned object give more information.
##'
##' @param a vector representing a `factor`, numeric, or alphabetically ordered character strings.  Censoring points have values of `-Inf`.
##' @param b like `a`.  If omitted, it copies `a`, representing nothing but uncensored values.  Censoring points have values of `Inf`.
##' @param precision when `a` and `b` are numeric, values may need to be rounded to avoid unpredictable behavior with \code{unique()} with floating-point numbers. Default is to 7 decimal places.
##' @param maxit maximum number of iterations allowed in the interval consolidation process when `cons='data'`
##' @param nponly set to `TRUE` to return a list containing the survival curve estimates before interval consolidation, using [icenReg::ic_np()]
##' @param cons set to `'none'` to not consolidate intervals when the survival estimate stays constant; this will likely cause a lot of trouble with zero cell probabilities during maximum likelihood estimation.  The default is to consolidate consecutive intervals.  Set `cons='data'` to change the raw data values to make observed intervals wider, in an iterative manner until no more consecutive tied survival estimates remain.
##' @param verbose set to `TRUE` to print information messages.  Set `verbose` to a number greater than 1 to get more information printed, such as the estimated survival curve at each stage of consolidation.
##' @return a 2-column integer matrix of class `"Ocens"` with an attribute `levels` (ordered), and if there are zero-width intervals arising from censoring, an attribute `upper` with the vector of upper limits.  Left-censored values are coded as `-Inf` in the first column of the returned matrix, and right-censored values as `Inf`.  When the original variables were `factor`s, these are factor levels, otherwise are numerically or alphabetically sorted distinct (over `a` and `b` combined) values.  When the variables are not factors and are numeric, other attributes `median`, `range`, and `npsurv` are also returned.  `median` is the median of the uncensored values on the origiinal scale.  `range` is a 2-vector range of original data values before adjustments to closest uncensored points.  `npsurv` is the estimated survival curve (with elements `time` and `surv`) from the `icenReg` package after any interval consolidation.  If the argument `npsurv=TRUE` was given, this `npsurv` list before consolidation is returned and no other calculations are done.  When the variables are factor or character, the median of the integer versions of variables for uncensored observations is returned as attribute `mid`.  A final attribute `freq` is the vector of frequencies of occurrences of all values.  `freq` aligns with `levels`.
##' @author Frank Harrell
##' @export
Ocens <- function(a, b=a, precision=7, maxit=10, nponly=FALSE,
                  cons=c('intervals', 'data', 'none'), verbose=FALSE) {
  cons <- match.arg(cons)
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

    # Since neither variable is a factor we can assume they are ordered
    # numerics.  Compute Turnbull intervals
    if(any(b < a)) stop('some values of b are less than corresponding a values')
    
    # Consider right-censored values to be in an open interval (a, Inf).
    # This causes the creation of a new category if beyond all uncensored obs.
    # Note that a and b are integers at this point
    j <- is.infinite(b)
    eps <- if(mul == 1e0) 0.001 else 0.1
    a[j] <- a[j] + eps
    # Similar for left-censored values
    j <- is.infinite(a)
    b[j] <- b[j] - eps

    ymed <- median(a[uncensored], na.rm=TRUE) * mul
    if(! requireNamespace('icenReg', quietly=TRUE)) stop('The icenReg package must be installed to use Ocens')
    fmi <- utils::getFromNamespace('findMaximalIntersections', 'icenReg')

    iter <- 0
    mto  <- function(x) diff(range(x)) > 0   # more than one distinct value
    repeat {
      iter <- iter + 1
      if(iter > maxit) stop('exceeded maxit=', maxit, ' iterations for pooling intervals')
      it  <- fmi(as.double(a), as.double(b))
      L   <- it$mi_l
      R   <- it$mi_r
      # The integer Y matrix produced by Ocens is the mappings of observations to the (L, R) Turnbull intervals
      # Indexes created by fmi start with 0, we bump them to 1
      ai  <- it$l_inds + 1L
      bi  <- it$r_inds + 1L
      if(verbose > 1) prn(mul * cbind(a, b, La=L[ai], Lb=L[bi], Ra=R[ai], Rb=R[bi]))
      dicen <- data.frame(a=a, b=b, grp=rep(1, length(a)))
      g <- icenReg::ic_np(cbind(a, b) ~ grp, data=dicen, B=c(1,1))  # bug prevents usage of matrix without formula
      # Note: icenReg::getSCurves() will not run
      s <- g$scurves[[1]]$S_curves$baseline
      k <- length(L) - 1L
      np <- list(time = L * mul, surv = s[1 : (k + 1)])
      if(nponly) return(np)
      if(length(np$time) != length(np$surv))
        warning('vector length mismatch in icenReg::ic_np result from npsurv=TRUE')
      if(verbose > 1) print(cbind(t=np$time, 'S(t)'=np$surv))
      if(cons == 'none') break
      s <- round(np$surv, 7)
      su <- unique(s[duplicated(s)])
      if(! length(su)) break

      if(cons == 'intervals') {
        # Consolidate intervals and interval definitions
        # Merge intervals having the same survival estimates
        # Unlike cons = 'data' this doens't change the data
        # First cluster row numbers in Turnbull interval table by s
        if(length(L) != length(s)) stop('program logic error in Ocens')
        # Construct consolidated intervals by mapping all intervals in a
        # cluster to the sequential cluster number
        us <- sort(unique(s), decreasing=TRUE)
        nt <- length(us)
        # Compute mapping of old table rows to new rows
        # old : 1 : length(L)
        new <- c(1, 1 + cumsum(diff(s) < 0))
        Ln <- Rn <- numeric(nt)
        for(i in 1 : nt) {
          s1 <- us[i]
          # Build new table Ln, Rn, knowing that s goes along with L
          j <- which(s == s1)
          Ln[i] <- min(L[j])
          Rn[i] <- max(R[j])
        }
        if(verbose) {
          cat('\nIntervals before consolidation\n\n')
          print(mul * cbind(L, R))
          cat('\nIntervals after consolidation\n\n')
          print(mul * cbind(Ln, Rn))
        }
        L <- Ln
        R <- Rn
        # Transform row numbers in raw data
        ai <- new[ai]
        bi <- new[bi]
        j  <- ! duplicated(s)
        np <- list(time = np$time[j], surv = np$surv[j])
        break
      }

      # Some consecutive intervals had the same information
      # For these code all the raw data as [lower, upper] where lower is the
      # minimum lower limit in the overlapping intervals, upper is the maximum upper limit
      # Compute distinct values of s that have > 1 Turnbull interval with that s value
      # Find original data corresponding to each su
      # Lookup s for each row of data
      S <- s[ai]

      for(ans in su) {
        j <- which(S == ans)
        if(! length(j)) stop('program logic error in Ocens')
        if(verbose) {
          cat('\nIntervals consolidated to give unique contributions to survival estimates and likelihood\n\nBefore:\n\n')
          print(cbind(a=mul * a[j], b=mul * b[j]))
        }
        aj <- a[j]
        bj <- b[j]
        l <- is.infinite(aj)
        r <- is.infinite(bj)
        ic <- (! l) & (! r) & (bj > aj)
        # Try only one remedy per group, using else if ...
        if(any(r) && ! any(l)) a[j[! l]] <- min(a[j[! l]])
        else if(any(r)) a[j[r]] <- min(aj)
        else if(any(l) && all(bj[l] == max(bj[! r]))) b[j[l]] <- min(bj[! r])
        else if(any(l)) b[r[l]] <- max(bj)
        else if((sum(ic) > 1) && (mto(a[j[ic]]) || mto(b[j[ic]]))) {
          a[j[ic]] <- min(aj[! l])
          b[j[ic]] <- max(bj[! r])
        }
        else if(any(ic)) {
          a[j] <- min(a[j])
          b[j] <- max(b[j])
        }
        if(verbose) {
          cat('\nAfter:\n\n')
          print(cbind(a=mul * a[j], b=mul * b[j]))
        }
      }
    }

    # freq is the count of number of observations
    # freq <- tabulate(ai[uncensored], nbins=length(u))
    freq <- tabulate(ai, nbins=max(ai))
    ai[is.infinite(a)] <- -Inf
    bi[is.infinite(b)] <-  Inf
    y    <- cbind(ai, bi)
    dimnames(y) <- list(NULL, NULL)

    return(structure(y,
                     class  = 'Ocens',
                     levels = L * mul,
                     upper  = if(any(L != R)) R * mul,
                     freq   = freq,
                     median = ymed,
                     range  = range(z),
                     npsurv = np))
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
