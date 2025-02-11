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

  f <- function(times=numeric(0), X=numeric(0),
                intercepts=numeric(0), slopes=numeric(0),
                info=numeric(0), values=numeric(0),
                interceptRef=integer(0), cumprob=cumprob, yname=NULL,
                conf.int=0)
  {
    np   <- length(slopes)
    k    <- length(intercepts)

    if(length(X) && (length(lp) > 1 || any(lp != 0))) stop('may not specify lp when X is given')
    if(conf.int && np > 0 && ! length(X))             stop('X must be specified when conf.int is used')

    # Compute linear predictor with no intercept used after saving original lp
    olp <- lp
    lp  <- if(np == 0) 0 else if(length(X)) matxv(X, slopes) else lp - intercepts[interceptRef]
    # Compute # distinct lps not counting possible varying intercepts
    ulp  <- length(unique(lp))
    cump <- eval(cumprob)    # cumprob is an expression defining a function
    lt   <- length(times)
    fvalues <- values[- length(values)]
    if(! lt) times <- fvalues

    # ints = indexes of intercepts corresponding to requested times, noting that the survival
    # function has estimates carried forward until the next distinct time
    # approx() will result in NA elements in ints for requested times outside range of values
    xrow  <- 1 : ulp
    w     <- expand.grid(time = times, Xrow=xrow)
    lps   <- lp[w$Xrow]
    ints  <- approx(fvalues, 1 : k, xout=w$time, method='constant')$y
    ki    <- length(ints)
    lpsi  <- intercepts[ints] + lps
    surv  <- cump(lpsi)
    xrow  <- w$Xrow
    w     <- data.frame(time=w$time, surv, Xrow=xrow)
    if(ulp == 1 || ! length(X)) w$Xrow <- NULL
    if(ulp > 1  && ! length(X)) w$lp   <- olp[xrow]
    row.names(w) <- NULL
    if(! conf.int) {
      r <- if(nrow(w) == 1) structure(w$surv, names=paste0('P(T>', w$time, ')'))
        else if(lt == 1)      w$surv
        else w
      return(r)
    }
    zcrit <- qnorm((1 + conf.int) / 2)
    
    # If there are Xs compute the column numbers in the information matrix that
    # correspond to just the columns of X being used
    # Index from 1 where the intercepts start

    idx <- integer(0)
    if(np > 0) {
      idx <- which(names(c(intercepts, slopes)) %in% colnames(X))
      X <- as.matrix(cbind(1, X))
    }
  se <- rep(NA, ki)
  for(ir in 1 : ki) {
    i  <- ints[ir]     # intercept number in play
    if(is.na(i)) next  # times value requested is out of range
    ix <- xrow[ir]     # row number of X in play
    if(np == 0) v <- infoMxop(info, i=i)
    else {
      Xi <- X[ix,, drop=FALSE]
      v  <- Matrix::diag(Xi %*% infoMxop(info, i=c(i, idx), B=t(Xi)))
    }
    se[ir] <- sqrt(v)
  }
  w$lower <- cump(lpsi - zcrit * se)
  w$upper <- cump(lpsi + zcrit * se)
  row.names(w) <- NULL
  w
  }

  cumprob <- object$famfunctions[1]

  formals(f) <-
    list(times=NULL, lp=0, X=numeric(0), 
                     intercepts=object$coef[1:ns],
                     slopes=object$coef[-(1 : ns)],
                     info=object$info.matrix, values=vals,
                     interceptRef=object$interceptRef,
                     cumprob=cumprob, yname=all.vars(object$terms)[1],
                     conf.int=0)
  f
}
