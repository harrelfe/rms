val.surv <- function(fit, newdata, S, est.surv, 
                     method=c('hare', 'smoothkm'),
                     censor, u, fun, lim, evaluate=100, pred, maxdim=5, 
                     ...)
{
  method  <- match.arg(method)
  isorm   <- ! missing(fit) && inherits(fit, 'orm')
  if(! missing(u)) {
    if(missing(fun)) {
      if(missing(fit))
        stop('must specify fit if u is specified and fun is not')
      fun <- if(inherits(fit, 'coxph')) function(p) log(-log(p))
        else if(isorm) eval(fit$famfunctions[2])
        else if(length(it <- survreg.distributions[[fit$dist]]$quantile)) it
        else if(length(it <- survreg.auxinfo[[fit$dist]]$quantile)) it
        else stop('fun is not given and it cannot be inferred')
    }
  }
  if(missing(S)) {
    S           <- fit[['y']]
    units       <- fit$units
    if(! length(S)) stop('when S is omitted you must use y=TRUE in the fit')
  }
  else units <- attr(S, 'units')

  if(NCOL(S) == 1) S <- Surv(S)
  if(! inherits(S, c('Surv', 'Ocens'))) stop('S must be a Surv or Ocens object')
  if(inherits(S, 'Ocens')) S <- Ocens2Surv(S)
  if(ncol(S) != 2) stop('S must be a right-censored Surv object')
  
  if(missing(est.surv))
    est.surv <- if(! missing(u)) {
      if(missing(newdata))
        survest(fit, times=u, conf.int=0)$surv
      else
        survest(fit, newdata, times=u, conf.int=0)$surv
    }
    else {
      if(missing(newdata))
        survest(fit, times=S[,1], what='parallel', conf.int=0)
      else
        survest(fit, newdata, times=S[,1], what='parallel', conf.int=0)
    }
  est.surv <- unname(est.surv)
  if(! missing(u)) {
    i <- ! is.na(est.surv + S[,1] + S[,2])
    est.surv   <- est.surv[i]
    S          <- S[i,, drop=FALSE]
    curtail    <- function(x) pmin(.9999, pmax(x, .0001))
    func       <- function(x) fun(curtail(x))
    xx         <- func(est.surv)

    f <- switch(method,
                hare     = polspline::hare(S[,1], S[,2], xx,
                                           maxdim=maxdim, ...),
                smoothkm = movStats(S ~ xx, times=u, melt=TRUE, ...) )
    if(missing(pred)) {
      if(missing(lim))
        lim <-
          datadist(est.surv)$limits[c('Low:prediction', 'High:prediction'),]
      pseq <- seq(lim[1], lim[2], length=evaluate)
    }
    else pseq <- pred
    
    if(method == 'hare') {
      actual    <- 1 - polspline::phare(u, func(est.surv), f)
      actualseq <- 1 - polspline::phare(u, func(pseq),     f)
    } else {
      actual    <- approx(f$xx, 1 - f$incidence, xout=func(est.surv))$y
      actualseq <- approx(f$xx, 1 - f$incidence, xout=func(pseq))$y
    }
    w <- structure(list(p=est.surv, actual=actual,
                        pseq=pseq, actualseq=actualseq,
                        u=u, fun=fun, n=nrow(S), d=sum(S[,2]),
                        units=units), class='val.survh')
    if(method == 'hare') w$harefit <- f else w$movstats <- f
    return(w)
  }
  
  n <- nrow(S)
  nac <- if(! missing(fit)) fit$na.action
  if(! missing(censor) && length(censor) > 1 && ! missing(fit)) {
    if(length(censor) > n && length(nac)) {
      ## Missing observations were deleted during fit
      j <- ! is.na(naresid(nac, censor))
      censor <- censor[j]
    }
    if(length(censor) != n)
      stop("length of censor does not match # rows used in fit")
  }
  
  est.surv.censor <- lp <- NULL
  if(! missing(censor)) {
    if(missing(fit))
      stop('fit must be specified when censor is specified')
    est.surv.censor <- if(missing(newdata))
      survest(fit, times=censor, what='parallel', conf.int=0)
    else
      survest(fit, newdata, times=censor, what='parallel', conf.int=0)
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
  if(length(x$harefit)) {
    cat('hare fit:\n\n')
    print(x$harefit)
  }
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
    uni <- Punits(x$units, adds=FALSE, default='unit')
    if(x$u != 1) uni <- paste0(uni, 's')
    lab <- paste('Probability of Surviving ', format(x$u), ' ', uni,
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
  S               <- x$S
  est.surv        <- x$est.surv
  censor.est.surv <- x$censor.est.surv
  what            <- match.arg(what)
  type            <- match.arg(type)
  
  n   <- length(est.surv)
  nac <- x$na.action
  
  if(! missing(group)) {
      if(length(group) > n && length(nac)) {
        ## Missing observations were deleted during fit
        j <- ! is.na(naresid(nac, est.surv))
        group <- group[j]
      }
      if(length(group) != n)
        stop("length of group does not match # rows used in fit")
      if(! is.factor(group)) group <- 
        if(is.logical(group) || is.character(group)) 
          as.factor(group) else cut2(group, g=g.group)
    }
  
  if(length(censor.est.surv)) {
    if(missing(group)) group <- rep(1, length(censor.est.surv))
    i <- S[, 2]==1
    group <- group[i]
    if(sum(i) < 2) stop('fewer than 2 uncensored observations')
    y <- switch(what,
                difference = 1-est.surv - .5*(1-censor.est.surv),
                ratio      = (1 - est.surv) / (.5 * (1 - censor.est.surv)))
    meanF <- tapply(1 - est.surv[i], group, mean, na.rm=TRUE)
    meanE <- tapply(.5 * (1 - censor.est.surv[i]), group, mean, na.rm=TRUE)
    res <- matrix(cbind(meanF, meanE), ncol=2,
                  dimnames=list(levels(group),
                    c('Mean F(T|T<C,X)', 'Expected')))
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
    nma <- ! is.na(est.surv + S[,2])
    est.surv <- est.surv[nma]
    S <- S[nma,, drop=FALSE]
    f <- survfitKM(factor(rep(1, length(est.surv))),
                   Surv(1 - est.surv,S[,2]),
                   se.fit = FALSE, conf.type = "none")
    tt <- c(0, f$time)
    ss <- c(1, f$surv)
    g <- ggplot(mapping=aes(x=tt, y=1 - ss)) + geom_step() +
      xlim(xlim) + ylim(ylim) + xlab(xlab) + ylab(ylab) +
      geom_abline(intercept = 0, slope = 1,
                  color = "black", linewidth=1, alpha=0.25)
    return(g)
  }
  
  nma <- ! (is.na(est.surv + S[,1] + S[,2]) | is.na(group))
  S        <- S[nma,, drop=FALSE]
  est.surv <- est.surv[nma]
  group    <- group[nma]
  
  ng <- length(lg <- levels(group))
  w <- NULL
  for(i in 1 : (ng + 1)) {
    s <- if(i==(ng + 1)) rep(TRUE,length(est.surv)) else group==lg[i]
    f <- survfitKM(factor(rep(1, sum(s))),
                   Surv(1 - est.surv[s], S[s, 2]),
                   se.fit = FALSE, conf.type = "none")
    w <- rbind(w, data.frame(group=c(lg, 'Overall')[i],
                             x=c(0, f$time), y=1 - c(1, f$surv)))
  }  
  g <- ggplot(w, aes(x, y, color=group)) + geom_step() +
        xlim(xlim) + ylim(ylim) + xlab(xlab) + ylab(ylab) +
        guides(color=guide_legend(title='')) +
        geom_abline(intercept = 0, slope = 1,
                    color = "black", linewidth=1, alpha=0.25)
  return(g)
}
