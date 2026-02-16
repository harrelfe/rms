##' Censored Ordinal Variable
##'
##' Combines two variables `a, b` into a 2-column matrix, preserving `label` and `units` attributes and converting character or factor variables into integers and added a `levels` attribute.  This is used to combine censoring points with regular points.  If both variables are already factors, their levels are distinctly combined starting with the levels for `a`.  Character variables are converted to factors.
##'
##' Left censored values will have `-Inf` for `a` and right-censored values will have `Inf` for `b`.  Interval-censored observations will have `b` > `a` and both finite.  For factor or character variables it only makes sense to have interval censoring.
##'
##' If there is no censoring, `a` is returned as an ordinary vector, with `label` and `units` attributes.
##'
##' @param a variable for first column
##' @param b variable for second column
##' @return a numeric matrix of class `Ocens`
##' @md
##' @author Frank Harrell
##' @export
Ocens <- function(a, b=a) {
  aname <- deparse(substitute(a))
  bname <- deparse(substitute(b))
  # If the arguments to Ocens were valid R names, use them
  name  <- if(aname == make.names(aname)) aname else if(bname == make.names(bname)) bname else ''

  i     <- ! is.na(a)
  if(any((! is.na(b)) != i)) stop('a and b must be NA on the same observations')

  ia    <- if(is.numeric(a)) is.finite(a) else ! is.na(a) # is.finite counts NAs also as FALSE
  ib    <- if(is.numeric(b)) is.finite(b) else ! is.na(b)
  uni   <- units(a)
  if(! length(uni) || uni == '') uni <- units(b)
  if(! length(uni)) uni <- ''
  lab   <- label(a)
  if(! length(lab) || lab == '') lab <- label(b)
  if(! length(lab)) lab <- ''

  if(is.character(a) + is.character(b) == 1) stop('neither or both of a and b should be character')
  if(is.factor(a) + is.factor(b) == 1)       stop('neither or both of a and b should be factor')
  if(is.numeric(a) + is.numeric(b) == 1)     stop('neither or both of a and b should be numeric')

  if(all(a[i] == b[i])) return(structure(a, label=lab, units=uni))

  lev <- NULL
  if(! is.numeric(a)) {
    if(is.character(a)) {
      lev <- sort(unique(c(a[ia], b[ib])))
      a   <- as.integer(factor(a, lev, lev))
      b   <- as.integer(factor(b, lev, lev))
    }
    else {   # factors
      alev <- levels(a)
      blev <- levels(b)
      # Cannot just pool the levels because ordering would not be preserved
      if(length(alev) >= length(blev)) {
        master <- alev
        other  <- blev
      } else {
        master <- blev
        other  <- alev
      }
      if(any(other %nin% master))
        stop('a variable has a level not found in the other variable')
      a   <- match(as.character(a), master)
      b   <- match(as.character(b), master)
      lev <- master
    }
  }
  structure(cbind(a=a, b=b), levels=lev, name=name, label=lab, units=uni, class='Ocens')
}


##' Recode Censored Ordinal Variable
##'
##' Creates a 2-column integer matrix that handles left- right- and interval-censored ordinal or continuous values for use in [rmsb::blrm()] and [orm()].  A pair of values `[a, b]` represents an interval-censored value known to be in the interval `[a, b]` inclusive of `a` and `b`. Left censored values are coded as `(-Infinity, b)` and right-censored as `(a, Infinity)`, both of these intervals being open at the finite endpoints. Open left and right censoring intervals are created by adding a small increment (subtracting for left censoring) to `a` or `b`. When this occurs at the outer limits, new ordinal categories will be created by `orm` to capture the real and unique information in outer censored values.  For example if the highest uncensored value is 10 and there is a right-censored value in the data at 10, a new category `10+` is created, separate from the category for `10`.  So it is assumed that if an exact value of 10 was observed, the pair of values for that observation would not be coded as `(10, Infinity)`.
##'
##' The intervals that drive the coding of the input data into numeric ordinal levels are the Turnbull intervals computed by the non-exported `findMaximalIntersections` function in the `icenReg` package, which handles all three types of censoring.  These are defined in the `levels` and `upper` attributes of the object returned by `Ocens`.  Sometimes consecutive Turnbull intervals contain the same statistical information likelihood function-wise, leading to the same survival estimates over two ore more consecutive intervals.  This leads to zero probabilities of involved ordinal values, preventing `orm` from computing a valid log-likeliihood.  A limited about of interval consolidation is done by `Ocens` to alleviate this problem.  Depending on the value of `cons` this consolidation is done by intervals (preferred) or by changing the raw data.  If `verbose=TRUE`, information about the actions taken is printed.
##'
##' When both input variables are `factor`s it is assumed that the one with the higher number of levels is the one that correctly specifies the order of levels, and that the other variable does not contain any additional levels.  If the variables are not `factor`s it is assumed their original values provide the orderings.  A left-censored point is is coded as having `-Inf` as a lower limit, and a right-censored point is coded as having `Inf` as an upper limit.   As with most censored-data methods, modeling functions assumes that censoring is independent of the response variable values that would have been measured had censoring not occurred.  `Ocens` creates a 2-column integer matrix suitable for ordinal regression.  Attributes of the returned object give more information.
##'
##' @param y an `Ocens` object, which is a 2-column numeric matrix, or a regular vector representing a `factor`, numeric, integer, or alphabetically ordered character strings.  Censoring points have values of `Inf` or `-Inf`.
##' @param precision when `y` columns are numeric, values may need to be rounded to avoid unpredictable behavior with \code{unique()} with floating-point numbers. Default is to 7 decimal places.  See [this](https://hbiostat.org/r/rms/unique-float/) for more details.
##' @param maxit maximum number of iterations allowed in the interval consolidation process when `cons='data'`
##' @param nponly set to `TRUE` to return a list containing the survival curve estimates before interval consolidation, using [icenReg::ic_np()]
##' @param cons set to `'none'` to not consolidate intervals when the survival estimate stays constant; this will likely cause a lot of trouble with zero cell probabilities during maximum likelihood estimation.  The default is to consolidate consecutive intervals.  Set `cons='data'` to change the raw data values to make observed intervals wider, in an iterative manner until no more consecutive tied survival estimates remain.
##' @param verbose set to `TRUE` to print information messages.  Set `verbose` to a number greater than 1 to get more information printed, such as the estimated survival curve at each stage of consolidation.
##' @return a 2-column integer matrix of class `"Ocens"` with an attribute `levels` (ordered), and if there are zero-width intervals arising from censoring, an attribute `upper` with the vector of upper limits.  Left-censored values are coded as `-Inf` in the first column of the returned matrix, and right-censored values as `Inf`.  When the original variables were `factor`s, these are factor levels, otherwise are numerically or alphabetically sorted distinct (over `a` and `b` combined) values.  When the variables are not factors and are numeric, other attributes `median`, `range`, `label`, and `npsurv` are also returned.  `median` is the median of the uncensored values on the origiinal scale.  `ranges` is a 3-element list, each element a 2-vector range.  The element named `y` is the range of original data values before adjustments.  The `u` element is a 2-vector range of uncensored values before adjustment, and the `c` element contains the lowest left censoring point and highest right-censored point.  Getting back to the main returned variables, `label` is the `label` attribute from the first of `a, b` having a label.  `npsurv` is the estimated survival curve (with elements `time` and `surv`) from the `icenReg` package after any interval consolidation.  If the argument `npsurv=TRUE` was given, this `npsurv` list before consolidation is returned and no other calculations are done.  When the variables are factor or character, the median of the integer versions of variables for uncensored observations is returned as attribute `mid`.  A final attribute `freq` is the vector of frequencies of occurrences of all values.  `freq` aligns with `levels`.  A `units` attribute is also included.  Finally there are two 3-vectors `Ncens1` and `Ncens2`, the first containing the original number of left, right, and interval-censored observations and the second containing the frequencies after altering some of the data.  For example, observations that are right-censored at or beyond the highest uncensored value are coded as uncensored to get the correct likelihood component in `orm.fit`.  When only right censoring is present and there are censored observations at or beyond the highest uncensored point, another attribute `rt_cens_beyond` is included in the returned list.  It has elements `newlevel` which is the numeric uncensored value assigned to these observations, and `range` which is a 2-vector containing the lowest and highest censored values beyond the last uncensored value.
##'
##' @author Frank Harrell
##' @export
Ocens2ord <- function(y, precision=7, maxit=10, nponly=FALSE,
                      cons=c('intervals', 'data', 'none'), verbose=FALSE) {
  cons <- match.arg(cons)
  # if(! inherits(y, 'Ocens')) stop('y must be an Ocens object')
  at <- attributes(y)

  if(NCOL(y) == 1) {
    a <- unclass(y)
    b <- a
  } else {
    a <- unclass(y)[, 1]
    b <- unclass(y)[, 2]
  }

  uni    <- at$units
  ylabel <- at$label
  if(ylabel == '') ylabel <- at$aname

  notna <- which(! is.na(a) & ! is.na(b))
  n <- length(a)
  A <- rep(NA_integer_, n)
  B <- rep(NA_integer_, n)
  if(length(notna) < length(a)) {
    a <- a[notna]
    b <- b[notna]
  }

  if(! length(at$levels)) {
    mul    <- 1e0
    z      <- c(a, b);  z <- z[is.finite(z)]
    yrange <- range(z)
    if(any(z %% 1 != 0)) {   # see recode2integer
      a   <- round(a * 10^precision)
      b   <- round(b * 10^precision)
      mul <- 10^-precision
    }
    uncensored <- a == b
    lc     <- is.infinite(a)
    rc     <- is.infinite(b)
    if(! any(uncensored)) stop('no uncensored observations')

    if(any(lc & rc)) stop('an observation has infinite values for both values')

    # Since neither variable is a factor we can assume they are ordered
    # numerics.  Compute Turnbull intervals
    if(any(b < a)) stop('some values of b are less than corresponding a values')

    urange <- range(a[uncensored]) * mul
    crange <- c(NA, NA)
    if(any(lc)) crange[1] <- min(b[lc]) * mul
    if(any(rc)) crange[2] <- max(a[rc]) * mul

    ymed   <- median(a[uncensored]) * mul

    # Compute original number of left, right, and interval-censored values
    ncen <- if(all(uncensored)) c(left=0, right=0, interval=0)
              else c(left=sum(lc), right=sum(rc),
                     interval=sum(is.finite(a) & is.finite(b) & a < b))

    if(sum(ncen) == 0) {
      u    <- sort(unique(a))
      y    <- match(a, u)
      freq <- tabulate(y, nbins=length(u))
      A[notna] <- y
      return(structure(cbind(a=A, b=A),
                       class  = 'Ocens',
                       levels = u * mul,
                       freq   = freq,
                       median = ymed,
                       ranges = list(y=yrange, u=urange, c=crange),
                       label  = ylabel,
                       units  = uni       )  )
      }

    eps <- if(mul == 1e0) 0.001 else 0.1

    # If only censored obs are right-censored, make simple adjustments
    # and compute Kaplan-Meier estimates

    min.outer.censored <- rt_cens_beyond <- NULL

    if(ncen[1] + ncen[3] == 0) {    # right censoring only
      # For obs censored at or beyond the last uncensored point, make a new
      # uncensored category at the minimum of such points, + eps

      maxu <- max(a[uncensored])
      i <- which(is.infinite(b) & a >= maxu)
      if(length(i)) {
        min.outer.censored <- min(a[i]) + eps
        rng           <- range(a[i])
        a[i]          <- min.outer.censored
        b[i]          <- a[i]
        uncensored[i] <- TRUE
        rt_cens_beyond <- list(newlevel = min.outer.censored * mul,
                               range    = rng * mul)
      }
      #?? Consider all other right censoring points as defining open intervals
      # a[! uncensored] <- a[! uncensored] + eps

      s <- km.quick(Surv(a, is.finite(b)), interval='>=')
      if(nponly) return(list(time=mul * s$time, surv=s$surv) )
      # u    <- c(sort(unique(a[uncensored])), min.outer.censored)
      u <- sort(unique(a[uncensored]))
      if(length(s$time) != length(u) || ! all.equal(u, s$time))
        stop('program logic error in Ocens: km.quick mismatch')
      y    <- match(a, u)
      # y may be an censored obs that has identical time to an uncensored one
      # This is OK
      y2   <- ifelse(is.infinite(b), b, y)
      # All non-matches should be censored values
      nm <- which(is.na(y) &  uncensored)
      if(any(nm)) {
        prn(min.outer.censored * mul)
        prn(u * mul)
        prn(nm)
        prn(cbind(y, y2, is.na(y), ! uncensored)[nm, ])
        stop('Ocens program logic error on non-matches')
      }
      freq <- tabulate(y[! is.na(y)], nbins=length(u))
      # Set all censored values < min.outer.censored to next smaller uncensored
      # values, and leave them censored
      # This is for censored values that were not tied with uncensored ones
      j   <- which(is.na(y))
      nl <- 0
      n  <- length(u)
      for(i in j) {
        below  <- which(u < a[i])
        if(length(below)) y[i] <- max(below)
        else {
          y[i] <- NA
          nl   <- nl + 1
        }
      }
      if(nl > 0) message(nl, ' observations are right-censored before any uncensored points.\n',
                          'These are set to NA.')

       s$time <- mul * s$time

      ncen2 <- if(all(uncensored)) c(left=0, right=0, interval=0)
              else c(left=sum(is.infinite(y)), right=sum(is.infinite(y2)),
                     interval=sum(is.finite(y) & is.finite(y2) & y < y2))


       A[notna] <- y
       B[notna] <- y2
       return(structure(cbind(a=A, b=B),
                        class   = 'Ocens',
                        levels  = u * mul,
                        freq    = freq,
                        median  = ymed,
                        ranges  = list(y=yrange, u=urange, c=crange),
                        rt_cens_beyond = rt_cens_beyond,
                        label   = ylabel,
                        units   = uni,
                        Ncens1  = ncen,
                        Ncens2  = ncen2,
                        npsurv  = s)          )
    }

    # What remains is left and interval censoring
    # Consider right-censored values to be in an open interval (a, Inf).
    # This causes the creation of a new category if beyond all uncensored obs.
    # Note that a and b are integers at this point
    j <- is.infinite(b)
    a[j] <- a[j] + eps
    # Similar for left-censored values
    j <- is.infinite(a)
    b[j] <- b[j] - eps

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
      s <- 1e-7 * round(np$surv * 1e7)
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

    ncen2 <- if(all(uncensored)) c(left=0, right=0, interval=0)
              else c(left=sum(is.infinite(ai)), right=sum(is.infinite(bi)),
                     interval=sum(is.finite(ai) & is.finite(bi) & (ai < bi)))

    A[notna] <- ai
    B[notna] <- bi
    y        <- cbind(a=A, b=B)
    dimnames(y) <- list(NULL, NULL)

    return(structure(y,
                     class   = 'Ocens',
                     levels  = L * mul,
                     upper   = if(any(L != R)) R * mul,
                     freq    = freq,
                     median  = ymed,
                     ranges  = list(y=yrange, u=urange, c=crange),
                     label   = ylabel,
                     units   = uni,
                     Ncens1  = ncen,
                     Ncens2  = ncen2,
                     npsurv  = np))
  }

  # Categorical variables as integers
  uncensored <- a == b
  if(! any(uncensored)) stop('no uncensored observations')
  if(any(b < a)) stop('some values of b are less than corresponding a values')
  freq <- tabulate(a[uncensored], nbins=length(at$levels))
  mid  <- quantile(a[uncensored], probs=.5, type=1L)
  A[notna] <- a
  B[notna] <- b
  # Categorical variables cannot be infinite, so no left or rt censoring
  ncen <- c(left=0, right=0, interval=sum(! uncensored))
  structure(cbind(a=A, b=B),
            class='Ocens', levels=at$levels, freq=freq, mid=mid,
            label=ylabel,
            units=uni,
            Ncens1=ncen, Ncens2=ncen)
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
##' @method as.data.frame Ocens
##' @export
as.data.frame.Ocens <- function(x, row.names = NULL, optional = FALSE, ...) {
  deb <- Fdebug('rmsdebug')
  nrows <- NROW(x)
  deb(nrows)
  row.names <- if(optional) character(nrows) else as.character(1:nrows)
  value <- list(x)
  deb(dim(value[[1]]))
  if(! optional) names(value) <- deparse(substitute(x))[[1]]
  deb(dim(value[[1]]))
  structure(value, row.names=row.names, class='data.frame')
}

##' Subset Method for `Ocens` Objects
##'
##' Subsets an `Ocens` object, preserving its special attributes.  Attributes are not updated.  In the future such updating should be implemented.
##' @title Ocens
##' @param x an `Ocens` object
##' @param ... the usual rows and columns specifiers
##' @param drop set to `FALSE` to not drop unneeded dimensions
##' @return new `Ocens` object or by default an unclassed vector if only one column of `x` is being kept
##' @author Frank Harrell
##' @md
##' @method [ Ocens
##' @export
'[.Ocens' <- function(x, ..., drop) {
  d  <- dim(x)
  at <- attributes(x)
  n  <- intersect(names(at), c('name', 'label', 'units', 'levels'))
  x  <- unclass(x)
  x  <- x[..., drop=FALSE]
  if(missing(drop)) drop <- NCOL(x) == 1
  if(drop) x <- drop(x)
  attributes(x) <- c(attributes(x), at[n])
  if(NCOL(x) == 2) class(x) <- 'Ocens'
  x
  }

##' is.na Method for Ocens Objects
##'
##' @param x an object created by `Ocens`
##'
##' @returns a logical vector whose length is the number of rows in `x`, with `TRUE` designating observations having one or both columns of `x` equal to `NA`
##' @method is.na Ocens
##' @export
##'
##' @md
##'
##' @examples
##' Y <- Ocens(c(1, 2, NA, 4))
##' Y
##' is.na(Y)
is.na.Ocens <- function(x) as.vector(rowSums(is.na(unclass(x))) > 0)

#' Ocens2Surv
#'
#' Converts an `Ocens` object to the simplest `Surv` object that works for the types of censoring that are present in the data.
#'
#' @param Y an `Ocens` object
#'
#' @returns a `Surv` object
#' @export
#' @md
#'
#' @examples
#' Y <- Ocens(1:3, c(1, Inf, 3))
#' Ocens2Surv(Y)
Ocens2Surv <- function(Y) {
  y  <- Y[, 1]
  y2 <- Y[, 2]

  su <- survival::Surv
  if(all(y == y2)) return(su(y))  # no censoring
  i <- which(is.finite(y) & is.finite(y2))
  w <- 1 * any(is.infinite(y)) + 2 * any(is.infinite(y2)) + 4 * any(y[i] != y2[i])
  if(w == 1)      su(y2, event=y == y2, type='left')
  else if(w == 2) su(y,  event=y == y2, type='right')
  else if(w == 4) su(y, event=rep(3, length(y)), time2=y2,       type='interval')
  else            su(y, time2=y2,       type='interval2')
}

##' print Method for Ocens Objects
##'
##' @param x an object created by `Ocens`
##' @param ivalues set to `TRUE` to print integer codes instead of character levels when original data were factors or character variables
##' @param digits number of digits to the right of the decimal place used in rounding original levels when `ivalues=FALSE`
##' @param ... ignored
##' @returns nothing
##' @method print Ocens
##' @export
##' @md
##'
##' @examples
##' Y <- Ocens(1:3, c(1, Inf, 3))
##' Y
##' print(Y, ivalues=TRUE)  # doesn't change anything since were numeric
print.Ocens <- function(x, ivalues=FALSE, digits=5, ...) {
  y   <- matrix(NA, nrow(x), ncol(x))   # to drop attributes of x
  y[] <- x
  a   <- y[, 1]
  b   <- y[, 2]
  nna <- ! is.na(a + b)
  ia  <- is.finite(a) & nna
  ib  <- is.finite(b) & nna
  ifa <- is.infinite(a) & nna
  ifb <- is.infinite(b) & nna
  lev <- attr(x, 'levels')
  if(! ivalues && length(lev) ) {
    a[ia] <- lev[a[ia]]
    b[ib] <- lev[b[ib]]
   }
  if(! length(lev)) {
    a <- round(a, digits)
    b <- round(b, digits)
  }
  intcens <- ia & ib & (b > a)
  a <- format(a)
  b <- format(b)
  z <- a
  z[ifa] <- paste0(b[ifa], '-')
  z[ifb] <- paste0(a[ifb], '+')
  z[intcens] <- paste0('[', a[intcens], ',', b[intcens], ']')
  print(z, quote=FALSE)
  invisible()
}

extractCodedOcens <- function(x, what=1, ivalues=FALSE, intcens=c('mid', 'low')) {
  intcens <- match.arg(intcens)
  lev   <- attr(x, 'levels')
  n <- nrow(x)
  a <- b <- integer(n)
  a[] <- x[, 1]   # gets rid of attributes
  b[] <- x[, 2]
  ia <- is.infinite(a)
  ib <- is.infinite(b)
  if(ivalues) {
    a <- a - 1
    b <- b - 1
  }
  else if(length(lev)) {
    a[! ia] <- lev[a[! ia]]
    b[! ib] <- lev[b[! ib]]
   }
  if(what == 2) return(cbind(a=a, b=b))

  ctype <- integer(n)
  ctype[ia]       <- 1                # left censoring
  ctype[ib]       <- 2                # right
  ctype[ctype == 0 & (a < b)] <- 3    # interval
  l <- ctype == 1
  r <- ctype == 2
  i <- ctype == 3

  y <- numeric(n)
  y[l] <- b[l]
  y[r] <- a[r]
  y[i] <- if(intcens == 'mid') 0.5 * (a + b)[i] else a[i]
  if(what == 1) return(y)
  list(a=a, b=b, y=y, ctype=ctype)
}

# Function determining TRUE/FALSE whether Y is known to be >= j
# a and b are results of Ocens2ord
# Returns NA if censoring prevents determining this
# Left censoring
#   Y >= j can be determined if b <= j
#   FALSE in this case
# Right censoring
#   Y >= j can be determined if a >= j
#   TRUE in this case
# Interval censoring
#   Y >= j can be determined if a >= j | b < j
#   TRUE if a >= j, FALSE if b < j
# Assumes that a and b run from 0 to k
geqOcens <- function(a, b, ctype, j) {
  z <- rep(NA, length(a))
  u <- ctype == 0
  l <- ctype == 1
  r <- ctype == 2
  i <- ctype == 3
  z[u] <- a[u] >= j
  z[l & b <= j] <- FALSE
  z[r & a >= j] <- TRUE
  z[i & a >= j] <- TRUE
  z[i & b <  j] <- FALSE
  z
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
