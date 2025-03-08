Survival.orm <- function(object, ...)
{
  ns   <- object$non.slopes
  vals <- as.numeric(object$yunique)
  up   <- object$yupper
  if(length(up)) {
      # When one limit is infinite use the other, otherwise take the midpoint
      vals <- ifelse(is.infinite(vals), up,
                     ifelse(is.infinite(up), vals, (vals + up) / 2.))
      vals <- as.numeric(vals)  # remove attributes
    }

  f <- function(times=numeric(0), lp=0, X=numeric(0),
                intercepts=numeric(0), slopes=numeric(0),
                info=numeric(0), values=numeric(0),
                interceptRef=integer(0), cumprob=cumprob, yname=NULL,
                yrange, conf.int=0, parallel=FALSE, forcedf=FALSE, zero=FALSE)
  {
    np   <- length(slopes)
    k    <- length(intercepts)

    if(length(X) && (length(lp) > 1 || any(lp != 0))) stop('may not specify lp when X is given')
    if(conf.int && np > 0 && ! length(X))             stop('X must be specified when conf.int is used')
    if(parallel && conf.int > 0)                      stop('conf.int must be 0 if parallel=TRUE')
    if(parallel && (length(lp) != length(times)))     stop('lp and times must have same length for parallel=TRUE')

    cump <- eval(cumprob)    # cumprob is an expression defining a function
    # Compute linear predictor with no intercept used after saving original lp
    olp <- lp
    lp  <- if(np == 0) 0 else if(length(X)) matxv(X, slopes) else lp - intercepts[interceptRef]

    # Compute # distinct lps not counting possible varying intercepts
    # TODO: watch out for different rows of X leading to same lp
    # No let's keep lp as-is
    i    <- TRUE    # ! duplicated(lp)
    olp  <- olp[i]
    ulp  <- lp[i]
    lulp <- length(ulp)
    lt   <- length(times)
    # Ordinal model is interested in P(Y >= highest value) but for survival
    # curves we want P(Y > last value) = 0
    fvalues <- values[- length(values)]
    zero <- zero && ! lt && (fvalues[1] != 0) && ! parallel
    # Get rid of any duplicate times specified by the user
    if(lt) {
      if(! parallel) times <- sort(unique(times))
    } else times <- fvalues

    # ints = indexes of intercepts corresponding to requested times, noting that the survival
    # function has estimates carried forward until the next distinct time
    # approx() will result in NA elements in ints for requested times outside range of values
    
    if(parallel) {
      ints <- approx(fvalues, 1 : k, xout=times, method='constant')$y
      r <- cump(intercepts[ints] + lp)
      attr(r, 'intercepts') <- NULL
      if(forcedf) r <- data.frame(time=times, surv=r, lp)
      return(r)
    }

    xrow  <- 1 : lulp
    w     <- expand.grid(time = times, Xrow=xrow)
    ints  <- approx(fvalues, 1 : k, xout=w$time, method='constant')$y
    j     <- which(is.na(ints))
    if(length(j)) {
      if(! lt) stop('program logic error: unexpected undefined intercepts')
      warning('Some times= requested are out of range of the data and are ignored:\n',
              paste(w$time[j], collapse=', '))
      w <- subset(w, time %nin% w$time[j])
      ints <- ints[! is.na(ints)]
    }

    # ints  <- pmin(ints, k)
    ki    <- length(ints)

    lps   <- ulp[w$Xrow]
    lpsi  <- intercepts[ints] + lps
    surv  <- cump(lpsi)
    xrow  <- w$Xrow
    w     <- data.frame(time=w$time, surv, Xrow=xrow)
    if(lulp == 1 || ! length(X)) w$Xrow <- NULL
    if(lulp > 1  && ! length(X)) w$lp   <- olp[xrow]
    row.names(w) <- NULL
    if(! conf.int) {
      # If only one thing is varying, return a vector, otherwise a data.frame unless forcedf=TRUE
      row.names(w) <- NULL
      if(zero) {
        z <- subset(w, time==min(w$time))
        z$time <- 0
        z$surv <- 1
        w <- rbind(z, w)
      }
      r   <- if(forcedf)   w
        else if(lt == 1)   w$surv
        else if(lulp == 1) structure(w$surv, names=paste0('P(T>', w$time, ')'))
        else               w
      return(r)
    }
  zcrit <- qnorm((1 + conf.int) / 2)
  idx   <- if(np > 0) (k  + 1) : (k  + np)
  idx2  <- if(np > 0) (ki + 1) : (ki + np)
  X     <- cbind(1e0, X)
  se    <- rep(NA, ki)
  # Invert entire information matrix (for all intercepts needed)
  # if < 3000 intercepts requested.  This will be faster than ki partial inverses.
  if(ki < 3000) inv <- infoMxop(info, i=c(ints, idx))  # may be a subset of intercepts
  for(ir in 1 : ki) {
    i  <- ints[ir]     # intercept number in play
#    if(is.na(i)) next  # times value requested is out of range
    ix <- xrow[ir]     # row number of X in play
    if(np == 0) v <- if(ki < 3000) inv[ir, ir] else infoMxop(info, i=i)
    else {
      Xi <- X[ix, , drop=FALSE]
      if(ki < 3000) {
        invir <- inv[c(ir, idx2), c(ir, idx2), drop=FALSE]
        v <- Xi %*% invir %*% t(Xi)
      }
      else v  <- Xi %*% infoMxop(info, i=c(i, idx), B=t(Xi))
    }
    se[ir] <- sqrt(v)
  }
  w$lower <- cump(lpsi - zcrit * se)
  w$upper <- cump(lpsi + zcrit * se)
  if(zero) {
    z <- subset(w, time==min(w$time))
    z$time  <- 0
    z$surv  <- 1
    z$lower <- NA
    z$upper <- NA
    w <- rbind(z, w)
    }
  row.names(w) <- NULL
  w
  }

  cumprob <- object$famfunctions[1]

  formals(f) <-
    list(times=NULL, lp=0, X=numeric(0), 
                     intercepts   = unname(object$coef[1:ns]),
                     slopes       = object$coef[-(1 : ns)],
                     info         = object$info.matrix, values=vals,
                     interceptRef = object$interceptRef,
                     cumprob      = cumprob,
                     yname        = object$yname,
                     yrange       = object$yrange,
                     conf.int=0, parallel=FALSE, forcedf=FALSE, zero=FALSE)
  f
}

utils::globalVariables('time')
