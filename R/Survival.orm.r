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
                info=numeric(0), values=numeric(0), rt_cens_beyond=NULL,
                interceptRef=integer(0), cumprob=cumprob, yname=NULL,
                ranges,
                conf.int=0, parallel=FALSE, forcedf=FALSE, zero=FALSE)
  {
    np   <- length(slopes)
    # Make it so P(Y > highest Y) = 0
    # A fake highest intercept is used for prediction, but for the variance
    # the last real intercept is used (see variable ints_curtailed below)
    intercepts <- c(intercepts, -Inf)
    k          <- length(intercepts)
    len_t      <- length(times)
    len_lp     <- length(lp)

    if(length(X) && (len_lp > 1 || any(lp != 0))) stop('may not specify lp when X is given')
    if(conf.int && np > 0 && ! length(X))         stop('X must be specified when conf.int is used')
    if(parallel && conf.int > 0)                  stop('conf.int must be 0 if parallel=TRUE')
    if(parallel && (len_lp != length(times)))     stop('lp and times must have same length for parallel=TRUE')

    cump <- eval(cumprob)    # cumprob is an expression defining a function
    # Compute linear predictor with no intercept used after saving original lp
    olp    <- lp
    lp     <- if(np == 0) 0
      else if(length(X)) matxv(X, slopes)
      else lp - intercepts[interceptRef]
    if(length(X)) len_lp <- length(lp)

    # Initially the thought was to reduce lp down to distinct values, but two different
    # values of X can have the same lp but different variances, so this idea was abandoned

    # Ordinal model is interested in P(Y >= highest value) but for survival
    # curves we want P(Y > last value) = 0
    zero <- zero && ! len_t && (values[1] != 0) && ! parallel
    # Get rid of any duplicate times specified by the user
    if(len_t) {
      if(! parallel) times <- sort(unique(times))
    } else {
      times <- values
      # Remove extra point added to get the right likelihood (parameters) and
      # add the last right-censored point.  Extra point added when there is at
      # least one uncensored point at or beyond the highest censored point, and
      # we want to extend the survival curve out to end of follow-up (last censored pt)
      if(length(rt_cens_beyond)) {
        times <- setdiff(times, rt_cens_beyond$newlevel)
        # if(rt_cens_beyond$range[2] > max(times))
        #   times <- c(times, rt_cens_beyond$range[2])  # add max uncensored time
      }
    }

    # ints = indexes of intercepts corresponding to requested times, noting that the survival
    # function has estimates carried forward until the next distinct time
    # ?? approx() will result in NA elements in ints for requested times outside range of values

    if(parallel) {
      if(length(values) != k)
        stop('logic error in Survival.orm with parallel=TRUE  k=', k,
             ' length(values)=', length(values))
      ints <- approx(values, 1 : k, xout=times, method='constant', rule=2)$y
      r <- cump(intercepts[ints] + lp)
      attr(r, 'intercepts') <- NULL
      if(forcedf) r <- data.frame(time=times, surv=r, lp)
      return(r)
    }

    if(length(values) != k)
      stop('logic error in Survival.orm  k=', k,
           ' length(values)=', length(values))

    xrow  <- seq_len(len_lp)
    # For a data frame with all linear predictor and time combinations
    w     <- expand.grid(time = times, Xrow=xrow)
    ts    <- w$time
    lps   <- lp[w$Xrow ]
    olps  <- olp[w$Xrow]
    # In ints, intercept indexes can be duplicated, e.g., when times are close together
    # but also for different X
    ints <- approx(values, 1 : k, xout=ts, method='constant', rule=2)$y

    lpsi <- lps + intercepts[ints]
    surv <- cump(lpsi)
    lpsi_last <- lps + intercepts[k - 1]   # duplicated for all time within x but OK

    # When an ints is NA, lpsi will be NA.  When we know the survival curve is
    # 0 or 1 we can override many of those.  This happens when requested times are
    # (1) at or above the last uncensored time and that time
    #     has no censored times at or above it (S(t) = 0)
    # (2) at or below the first uncensored time and that time has no censored times
    #     at or below it (S(t) = 1)
    # Also if there are censored points beyond the last uncensored point, the survival
    # curve is carried forward from the last value getting
    # an intercept, to the outer censoring point

    min_uncens <- ranges$u[1]
    max_uncens <- ranges$u[2]
    surv[ts < min_uncens] <- 1.0
    LP <- ifelse(ts < min_uncens, Inf, lpsi)
    if(length(rt_cens_beyond)) {   # rt-censored points at or beyond highest uncensored pt
      rng <- rt_cens_beyond$range
      max_fu <- rng[2]
      surv[ts > max_fu] <- NA
      LP   <- ifelse(ts > rng[1] & ts <= max_fu, lpsi_last, LP)
      LP[ts > max_fu]   <- NA
      surv <- cump(LP)
    } else {
      surv[ts >= max_uncens] <- 0.0
    }

    xrow  <- w$Xrow
    w     <- data.frame(time=ts, surv, Xrow=xrow)
    if(len_lp == 1 || ! length(X)) w$Xrow <- NULL
    if(len_lp > 1  && ! length(X)) w$lp   <- olps
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
      r   <- if(forcedf)     w
        else if(len_t  == 1) w$surv
        else if(len_lp == 1) structure(w$surv, names=paste0('P(T>', w$time, ')'))
        else                 w
      return(r)
    }

  zcrit          <- qnorm((1 + conf.int) / 2)
  len_ints       <- length(ts)
  ints_curtailed <- pmin(ints, k - 1)

  # k is one plus number of intercepts to subtract 1 from it
  idx   <- if(np > 0) k : (k  + np - 1)
  idx2  <- if(np > 0) (len_ints + 1) : (len_ints + np)
  X     <- cbind(1e0, X)
  se    <- rep(NA, len_ints)

  # Invert entire information matrix (for all intercepts needed)
  # if < 3000 intercepts requested.  This will be faster than len_ints partial inverses.
  # I verified that infoMxop, when computing the inverse matrix for elements in i
  # properly handles duplicates in i

  if(len_ints < 3000) inv <- infoMxop(info, i=c(ints_curtailed, idx))
  # may be a subset of intercepts
  for(ir in 1 : len_ints) {
    i  <- ints_curtailed[ir]  # intercept number in play
    ix <- xrow[ir]            # row number of X in play
    if(np == 0) v <- if(len_ints < 3000) inv[ir, ir] else infoMxop(info, i=i)
    else {
      Xi <- X[ix, , drop=FALSE]
      if(len_ints < 3000) {
        invir <- inv[c(ir, idx2), c(ir, idx2), drop=FALSE]
        v <- Xi %*% invir %*% t(Xi)
      }
      else v  <- Xi %*% infoMxop(info, i=c(i, idx), B=t(Xi))
    }
    se[ir] <- sqrt(v)
  }
  w$lower <- cump(LP - zcrit * se)   # LP does carry forward
  w$upper <- cump(LP + zcrit * se)

  w$lower[w$surv == 1.0] <- 1.0   # matches survival::survfit convention
  w$upper[w$surv == 1.0] <- 1.0
  w$lower[is.na(w$surv) | w$surv == 0.0] <- NA
  w$upper[is.na(w$surv) | w$surv == 0.0] <- NA

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
                     intercepts     = unname(object$coef[1:ns]),
                     slopes         = object$coef[-(1 : ns)],
                     info           = object$info.matrix,
                     values         = vals,
                     rt_cens_beyond = object$rt_cens_beyond,
                     interceptRef   = object$interceptRef,
                     cumprob        = cumprob,
                     yname          = object$yname,
                     ranges         = object$ranges,
                     conf.int=0, parallel=FALSE, forcedf=FALSE, zero=FALSE)
  f
}

utils::globalVariables('time')
