val.surv <- function(fit, newdata, S, est.surv, censor,
                     u, fun, lim, evaluate=100, pred, maxdim=5, ...)
{
  usehare <- !missing(u)
  if(usehare) {
    require(polspline) || stop('must have polspline installed')
    if(missing(fun)) {
      if(missing(fit))
        stop('must specify fit if u is specified and fun is not')
      fun <- if(inherits(fit, 'coxph')) function(p) log(-log(p))
      else if(length(it <- survreg.distributions[[fit$dist]]$quantile)) it
      else if(length(it <- survreg.auxinfo[[fit$dist]]$quantile)) it
      else stop('fun is not given and it cannot be inferred')
    }
  }
  if(missing(S)) {
    S     <- fit$y
    units <- fit$units
    if(!length(S)) stop('when S is omitted you must use y=TRUE in the fit')
    if(!usehare) {
      itrans <- survreg.distributions[[fit$dist]]$itrans
      S[,1]  <- itrans(S[,1])
    }
  }
  else units <- attr(S, 'units')
  
  if(!any(attr(S,'class')=='Surv')) stop('S must be a Surv object')
  if(ncol(S) != 2) stop('S must be a right-censored Surv object')
  if(missing(est.surv))
    est.surv <- if(usehare) {
      if(missing(newdata))
        survest(fit, times=u)$surv
      else
        survest(fit, newdata, times=u)$surv
    }
    else {
      if(missing(newdata))
        survest(fit, times=S[,1], what='parallel')
      else
        survest(fit, newdata, times=S[,1], what='parallel')
    }
  if(usehare) {
    i <- !is.na(est.surv + S[,1] + S[,2])
    est.surv   <- est.surv[i]
    S <- S[i,]
    curtail <- function(x) pmin(.9999, pmax(x, .0001))
    f <- hare(S[,1], S[,2], fun(curtail(est.surv)), maxdim=maxdim, ...)
    if(missing(pred)) {
      if(missing(lim))
        lim <-
          datadist(est.surv)$limits[c('Low:prediction','High:prediction'),]
      pseq <- seq(lim[1], lim[2], length=evaluate)
    }
    else pseq <- pred
    
    actual    <- 1 - phare(u, fun(curtail(est.surv)), f)
    actualseq <- 1 - phare(u, fun(curtail(pseq)),     f)
    w <- structure(list(harefit=f, p=est.surv, actual=actual,
                        pseq=pseq, actualseq=actualseq,
                        u=u, fun=fun, n=nrow(S), d=sum(S[,2]),
                        units=units), class='val.survh')
    return(w)
  }
  
  n <- nrow(S)
  nac <- if(!missing(fit)) fit$na.action
  if(!missing(censor) && length(censor) > 1 && !missing(fit)) {
    if(length(censor) > n && length(nac)) {
      ## Missing observations were deleted during fit
      j <- !is.na(naresid(nac, censor))
      censor <- censor[j]
    }
    if(length(censor) != n)
      stop("length of censor does not match # rows used in fit")
  }
  
  est.surv.censor <- lp <- NULL
  if(!missing(censor)) {
    if(missing(fit))
      stop('fit must be specified when censor is specified')
    est.surv.censor <- if(missing(newdata))
      survest(fit, times=censor, what='parallel')
    else
      survest(fit, newdata, times=censor, what='parallel')
    if(mc <- sum(is.na(est.surv.censor)))
      warning(paste(mc, 'observations had missing survival estimates at censoring time'))
    lp <- if(missing(newdata)) predict(fit, type='lp')
    else
      predict(fit, newdata, type='lp')
  }
  
  if(length(est.surv) != n)
    stop('length of est.surv must equal number of rows in S')
  structure(list(S=S, est.surv=est.surv,
                 censor.est.surv=if(length(est.surv.censor))
                 est.surv.censor,
                 lp=if(length(lp))lp,
                 na.action=nac), class='val.surv')
}

print.val.survh <- function(x, ...) {
  cat('\nValidation of Predicted Survival at Time=', format(x$u),
      '\tn=', x$n, ', events=', x$d, '\n\n')
  cat('hare fit:\n\n')
  print(x$harefit)
  cat('\nFunction used to transform predictions:\n')
  cat(paste(format(x$fun), collapse=' '))
  error <- abs(x$p - x$actual)
  er <- c(mean(error, na.rm=TRUE), quantile(error, .9, na.rm=TRUE))
  erf <- format(round(er, 4))
  cat('\n\nMean absolute error in predicted probabilities:',
      erf[1],'\n')
  cat('0.9 Quantile of absolute errors               :', erf[2], '\n')
  er
}

plot.val.survh <- function(x, lim, xlab, ylab,
                           riskdist=TRUE, add=FALSE,
                           scat1d.opts=list(nhistSpike=200), ...)
  {
    if(missing(lim)) lim <- range(c(x$pseq, x$actualseq), na.rm=TRUE)
    units <- if(x$u == 1) x$units
    else
      paste(x$units, 's', sep='')
    lab <- paste('Probability of Surviving ', format(x$u), ' ', units,
                 sep='')
    if(add)
      lines(x$pseq, x$actualseq, ...)
    else
      plot(x$pseq, x$actualseq, type='l', xlim=lim, ylim=lim,
           xlab=if(missing(xlab)) paste('Predicted', lab) else xlab,
           ylab=if(missing(ylab)) paste('Actual',    lab) else ylab)
    abline(a=0, b=1, lty=2)
    if(riskdist) do.call('scat1d', c(list(x=x$p), scat1d.opts))
    invisible()
  }
    


plot.val.surv <- function(x, group, g.group=4,
                          what=c('difference','ratio'),
                          type=c('l','b','p'),
                          xlab, ylab, xlim, ylim,
                          datadensity=TRUE, ...)
{
  S <- x$S
  est.surv <- x$est.surv
  censor.est.surv <- x$censor.est.surv
  what <- match.arg(what)
  type <- match.arg(type)
  
  n <- length(est.surv)
  nac <- x$na.action
  
  if(!missing(group)) {
      if(length(group) > n && length(nac)) {
        ## Missing observations were deleted during fit
        j <- !is.na(naresid(nac, est.surv))
        group <- group[j]
      }
      if(length(group) != n)
        stop("length of group does not match # rows used in fit")
      if(!is.factor(group)) group <- 
        if(is.logical(group) || is.character(group)) 
          as.factor(group) else cut2(group, g=g.group)
    }
  
  if(length(censor.est.surv)) {
    if(missing(group)) group <- rep(1, length(censor.est.surv))
    i <- S[,2]==1
    group <- group[i]
    if(sum(i)<2) stop('fewer than 2 uncensored observations')
    y <- switch(what,
                difference=1-est.surv - .5*(1-censor.est.surv),
                ratio=(1 - est.surv) / (.5 * (1 - censor.est.surv)))
    meanF <- tapply(1 - est.surv[i], group, mean, na.rm=TRUE)
    meanE <- tapply(.5*(1-censor.est.surv[i]), group, mean, na.rm=TRUE)
    res <- matrix(cbind(meanF,meanE), ncol=2,
                  dimnames=list(levels(group),
                    c('Mean F(T|T<C,X)','Expected')))
    lp <- x$lp
    lp <- lp[i]; y <- y[i]
    if(missing(xlab)) xlab <- 'Linear Predictor'
    if(missing(ylab)) ylab <-
      switch(what,
             difference='F(T|X,T<=C) - .5F(C|X)',
             ratio='F(T|X,T<=C)/.5F(C|X), C=Censoring Time')
    if(missing(xlim)) xlim <- range(lp, na.rm=TRUE)
    if(missing(ylim)) ylim <- range(y,  na.rm=TRUE)
    if(type=='l') plsmo(lp, y, group=group, xlab=xlab, ylab=ylab,
         xlim=xlim, ylim=ylim,
         datadensity=datadensity, trim=0, ...)
    else {
      plot(lp, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
      if(type!='p') plsmo(lp, y, group=group, add=TRUE,
           xlab=xlab, ylab=ylab, datadensity=datadensity, trim=0, ...)
    }
    abline(h=if(what=='difference')0 else 1, lty=2)
    return(res)
  }
  
  if(missing(xlab)) xlab <- 'Predicted Pr[T <= observed T]'
  if(missing(ylab)) ylab <- 'Fraction <= x'
  if(missing(xlim)) xlim <- 0:1
  if(missing(ylim)) ylim <- 0:1
  
  if(missing(group)) {
    nma <- !is.na(est.surv + S[,2])
    est.surv <- est.surv[nma]
    S <- S[nma,,drop=FALSE]
    f <- survfitKM(factor(rep(1,length(est.surv))),
                   Surv(1 - est.surv,S[,2]),
                   se.fit = FALSE, conf.type = "none")
    tt <- c(0, f$time)
    ss <- c(1, f$surv)
    plot(tt, 1-ss, xlab=xlab, ylab=ylab, type='s',
         xlim=xlim, ylim=ylim, ...)
    abline(a=0, b=1, lty=2)
    return(invisible())
  }
  
  nma <- !(is.na(est.surv + S[,1] + S[,2]) | is.na(group))
  S <- S[nma,,drop=FALSE]
  est.surv <- est.surv[nma]
  
  ng <- length(lg <- levels(group))
  curves <- vector('list',ng+1)
  names(curves) <- c(lg, 'Overall')
  for(i in 1:(ng+1)) {
    s <- if(i==(ng+1)) rep(TRUE,length(est.surv)) else group==lg[i]
    f <- survfitKM(factor(rep(1,sum(s))),
                   Surv(1 - est.surv[s], S[s,2]),
                   se.fit = FALSE, conf.type = "none")
    curves[[i]] <- list(x=c(0,f$time), y=1-c(1,f$surv))
  }  
  labcurve(curves, pl=TRUE, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, ...)
  abline(a=0,b=1,lty=2)
  invisible()
}
