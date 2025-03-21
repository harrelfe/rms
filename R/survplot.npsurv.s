survplot.npsurv <-
  function(fit, xlim, 
           ylim, xlab, ylab, time.inc, state=NULL,
           conf=c("bands", "bars", "diffbands", "none"), mylim=NULL,
           add=FALSE, label.curves=TRUE, abbrev.label=FALSE, levels.only=FALSE,
           lty, lwd=par('lwd'),
           col=1, col.fill=gray(seq(.95, .75, length=5)),
           loglog=FALSE, fun, n.risk=FALSE, aehaz=FALSE, times=NULL,
           logt=FALSE, dots=FALSE, dotsize=.003, grid=NULL,
           srt.n.risk=0, sep.n.risk=.056, adj.n.risk=1,
           y.n.risk, cex.n.risk=.6, cex.xlab=par('cex.lab'), cex.ylab=cex.xlab,
           pr=FALSE, ...) {
    
    conf     <- match.arg(conf)
    polyg    <- ordGridFun(grid=FALSE)$polygon
    conf.int <- fit$conf.int
    if(!length(conf.int) | conf == "none") conf.int <- 0
    opar <- par(c('mar', 'xpd'))
    on.exit(par(opar))

    cylim <- function(ylim)
      if(length(mylim)) c(min(ylim[1], mylim[1]), max(ylim[2], mylim[2]))
    else ylim
    
    fit.orig <- fit
    units   <- Punits(fit$units, adds=FALSE, upfirst=TRUE, default='day')
    maxtime <- fit$maxtime
    if(! length(maxtime)) maxtime <- max(fit$time)
    mintime <- min(fit$time, 0)
    pret    <- pretty(c(mintime, maxtime))
    maxtime <- max(pret)
    mintime <- min(pret)
    if(missing(time.inc)) {
      time.inc <- switch(units, Day=30, Month=1, Year=1,
                                (maxtime - mintime) / 10)
      if(time.inc > maxtime) time.inc <- (maxtime - mintime) / 10
    }
    if(n.risk && ! length(fit$n.risk)) {
      n.risk <- FALSE
      warning("fit does not have number at risk\nIs probably from a parametric model\nn.risk set to FALSE")
    }

    mstate <- inherits(fit, 'survfitms')
    if(mstate) {
      ## Multi-state model for competing risks
      if(missing(fun)) fun <- function(y) 1 - y
      if(missing(state)) stop('state must be given when response is a multi-state/competing risk object from Surv()')
      if(length(state) != 1) stop('at present state can only be a single state')
      states <- fit$states
      if(state %nin% states) stop(paste('state is not in',
                                        paste(states, collapse=', ')))
    }
    
    trans <- loglog || mstate || ! missing(fun)
    if(missing(ylab))
      ylab <- 
       if(loglog) "log(-log Survival Probability)"
       else if(mstate) paste('Cumulative Incidence of', upFirst(state))
       else if(trans) ""
       else "Survival Probability"

    if(loglog) fun <- function(y) logb(-logb(ifelse(y == 0 | y == 1, NA, y)))
     else if(! trans) fun <- function(y) y
    
    if(missing(xlab))
      xlab <- if(logt) paste("log Follow-up Time in ", units, "s", sep="")
      else labelPlotmath('Follow-up Time', paste(fit$units, 's', sep=''))
      ## else labelPlotmath(fit$time.label, fit$units)
  
    if(missing(xlim)) 
      xlim <- if(logt) logb(c(maxtime / 100, maxtime)) else c(mintime, maxtime)

    convert <- if(mstate) {
      istate    <- match(state, states)
      conv <- function(f, istate) {
        f$surv    <- 1 - f$pstate [, istate]
        f$lower   <- 1 - f$lower  [, istate]
        f$upper   <- 1 - f$upper  [, istate]
        f$std.err <-     f$std.err[, istate]
        icens     <- which(states == '(s0)')
        if(! length(icens))
          stop('Program logic error: did not find (s0) column with competing risks')
        f$n.risk  <- f$n.risk[, icens]
        if(all(f$n.risk == 0))
          stop('program logic error: all n.risk are zero')
        f
      }
      formals(conv) <- list(f=NULL, istate=istate)
      conv } else function(f) f

  fit <- convert(fit)
  
  origsurv <- fit$surv
  if(trans) {
    fit$surv <- fun(fit$surv)
    fit$surv[is.infinite(fit$surv)] <- NA
    ##  handle e.g. logit(1) - Inf would mess up ylim in plot()
    if(conf.int > 0) {
      fit$lower <- fun(fit$lower)
      fit$upper <- fun(fit$upper)
      fit$lower[is.infinite(fit$lower)] <- NA
      fit$upper[is.infinite(fit$upper)] <- NA
      if(missing(ylim))
        ylim <- cylim(range(c(fit$lower, fit$upper), na.rm=TRUE))
    }
    else if(missing(ylim)) ylim <- cylim(range(fit$surv, na.rm=TRUE))
  }
  else if(missing(ylim)) ylim <- c(0, 1)

  if(length(grid)) {
    dots <- FALSE
    if(is.logical(grid)) grid <- if(grid) gray(.8) else NULL
  }
  if(logt | trans) { dots <- FALSE; grid <- NULL }

  olev <- slev <- names(fit$strata)
  if(levels.only) slev <- gsub('.*=', '', slev)
  sleva <- if(abbrev.label) abbreviate(slev) else slev
  ns <- length(slev)
  slevp <- ns > 0

  labelc <- is.list(label.curves) || label.curves
  if(!slevp) labelc <- FALSE

  ns <- max(ns, 1)

  y <- 1 : ns
  stemp <- if(ns == 1) rep(1, length(fit$time))
  else rep(1:ns, fit$strata)

  if(n.risk | (conf.int > 0 & conf == "bars")) {
    stime <- seq(mintime, maxtime, time.inc)
    v <- convert(summary(fit.orig, times=stime, print.it=FALSE))
    v$surv  <- fun(v$surv)
    v$lower <- fun(v$lower)
    v$upper <- fun(v$upper)
    vs <- if(ns > 1) as.character(v$strata)
    ## survival:::summary.survfit was not preserving order of strata levels
  }
  xd <- xlim[2] - xlim[1]
  yd <- ylim[2] - ylim[1]
  
  if(n.risk && !add) {
    mar <- opar$mar
    if(mar[4] < 4) {mar[4] <- mar[4] + 2; par(mar=mar)}
  }
  
  ## One curve for each value of y, excl style used for C.L.
  lty <- if(missing(lty)) seq(ns + 1)[-2] else rep(lty, length=ns)
  lwd <- rep(lwd, length=ns)
  col <- rep(col, length=ns)

  if(conf == 'diffbands' && ns < 2) conf <- 'bands'
  if(labelc || conf %in% c('bands', 'diffbands')) curves <- vector('list', ns)
  Tim <- Srv <- list()
  
  par(xpd=NA)

  nevents <- totaltime <- numeric(ns)
  cuminc  <- character(ns)
  for(i in 1 : ns) {
    st <- stemp == i
    time         <- fit$time[st]
    surv         <- fit$surv[st]
    osurv        <- origsurv[st]
    ## nevents[i]   <- sum(fit$n.event[st])
    ## nrsk         <- fit$n.risk[st]
    ## neachtime    <- c(- diff(nrsk), min(nrsk))
    ## totaltime[i] <- sum(neachtime * time)
    nevents[i] <- if(mstate) {
                    if(ns == 1) fit$numevents[, state]
                    else fit$numevents[olev[i], state]
    } else {
      if(ns == 1) fit$numevents else fit$numevents[olev[i]]
    }
    totaltime[i] <- if(ns == 1) fit$exposure else fit$exposure[olev[i]]
    if(length(times)) {
      cumi <- 1. - approx(time, osurv, xout=times, method='constant')$y
      noun <- units %in% c('', ' ')
      cuminc[i]   <- paste('Cum. inc.@ ',
                            if(noun) 't=',
                            paste(times, collapse=','),
                            if(! noun) paste(' ', units, sep=''),
                            ':', paste(round(cumi, 3), collapse=','),
                            sep='')
    }
      
    if(logt) time <- logb(time)
    s <- !is.na(time) & (time >= xlim[1])
    if(i==1 & !add) {
      plot(time, surv, xlab='', xlim=xlim,
           ylab='', ylim=ylim, type="n", axes=FALSE)
      mgp.axis(1, at=if(logt) pretty(xlim) else
               seq(xlim[1], max(pretty(xlim)), time.inc),
               labels=TRUE, axistitle=xlab, cex.lab=cex.xlab)
      mgp.axis(2, at=pretty(ylim), axistitle=ylab, cex.lab=cex.ylab)

      if(dots || length(grid)) {
        xlm <- pretty(xlim)
        xlm <-   c(xlm[1], xlm[length(xlm)])
        xp  <- seq(xlm[1], xlm[2], by=time.inc)
        yd  <- ylim[2] - ylim[1]
        if(yd <= .1) yi <- .01
        else if(yd <=.2 ) yi <- .025
        else if(yd <=.4 ) yi <- .05
        else yi <- .1
        yp <- seq(ylim[2],
                  ylim[1] + if(n.risk && missing(y.n.risk)) yi else 0,
                  by = -yi)
        if(dots)
          for(tt in xp)
            symbols(rep(tt, length(yp)), yp,
                    circles=rep(dotsize, length(yp)),
                    inches=dotsize, add=TRUE)
        else abline(h=yp, v=xp, col=grid, xpd=FALSE)
      }
    }
    tim <- time[s]; srv <- surv[s]
    if(conf.int > 0 && conf == "bands") {
      blower <- fit$lower[st][s]
      bupper <- fit$upper[st][s]
    }
    ##don't let step function go beyond x-axis -
    ##this cuts it off but allows step to proceed axis end
    if(max(tim) > xlim[2]) {
      srvl <- srv[tim <= xlim[2] + 1e-6]
      ## s.last <- min(srv[tim<=xlim[2]+1e-6]) #not work w/fun
      s.last <- srvl[length(srvl)]
      k <- tim < xlim[2]
      tim <- c(tim[k], xlim[2])
      srv <- c(srv[k], s.last)
      if(conf.int > 0 && conf == "bands") {
        low.last <- blower[time <= xlim[2] + 1e-6]
        low.last <- low.last[length(low.last)]
        up.last  <- bupper[time <= xlim[2] + 1e-6]
        up.last  <- up.last[length(up.last)]
        blower   <- c(blower[k], low.last)
        bupper   <- c(bupper[k], up.last)
      }
    }
    if(logt) {
      if(conf %nin% c('bands', 'diffbands'))
        lines(tim, srv, type="s", lty=lty[i], col=col[i], lwd=lwd[i])
      if(labelc || conf %in% c('bands', 'diffbands'))
        curves[[i]] <- list(tim, srv)
    }
	  else {
      xxx <- c(mintime, tim)
      yyy <- c(fun(1), srv)
      if(conf %nin% c('bands', 'diffbands'))
        lines(xxx, yyy, type="s", lty=lty[i], col=col[i], lwd=lwd[i])
      if(labelc || conf %in% c('bands', 'diffbands'))
        curves[[i]] <- list(xxx, yyy)
    }
    if(pr) {
      zest <- rbind(time[s], surv[s])
      dimnames(zest) <- list(c("Time", "Survival"),
                             rep("", sum(s)))
      if(slevp)cat("\nEstimates for ", slev[i], "\n\n")
      print(zest, digits=3)
    }
    if(conf.int > 0) {
      if(conf == 'bands') {
        if(logt)
          polyg(x = c(tim, max(tim), rev(tim)),
                y = c(blower, rev(bupper), max(bupper)),
                col =  col.fill[i], type='s')
        else
          polyg(x = c(max(tim), tim, rev(c(tim, max(tim)))),
                y = c(fun(1), blower, rev(c(fun(1), bupper))),
                col = col.fill[i], type = "s")
      }
      else if(conf == 'diffbands')
        survdiffplot(fit.orig, conf=conf, fun=fun, convert=convert, xlim=xlim)

      else {
        j <- if(ns == 1) TRUE else vs == olev[i]
        tt <- v$time[j]  #may not get predictions at all t
        ss <- v$surv[j]
        lower <- v$lower[j]
        upper <- v$upper[j]
        if(logt) tt <- logb(ifelse(tt == 0, NA, tt))
        tt <- tt + xd * (i - 1) * .01
        errbar(tt, ss, upper, lower, add=TRUE, lty=lty[i],
               col=col[i])
      }
    }
    
    if(n.risk) {
      j <- if(ns == 1) TRUE else vs == olev[i]
      tt <- v$time[j]
      nrisk <- v$n.risk[j]
      tt[1] <- xlim[1]  #was xd*.015, .030, .035
      if(missing(y.n.risk)) y.n.risk <- ylim[1]
      if(y.n.risk == 'auto') y.n.risk <- - diff(ylim) / 3
      yy <- y.n.risk + yd * (ns - i) * sep.n.risk
      nri <- nrisk
      nri[tt > xlim[2]] <- NA
      text(tt[1], yy, nri[1], cex=cex.n.risk,
           adj=adj.n.risk, srt=srt.n.risk)
      if (length(nri) > 1)
        text(tt[-1], yy, nri[-1], cex=cex.n.risk, adj=1)
      if(slevp) text(xlim[2] + xd * .025,
                    yy, adj=0, sleva[i], cex=cex.n.risk)
    }
  }
  
  if(conf %in% c('bands', 'diffbands'))
    for(i in 1:ns)
      lines(curves[[i]][[1]], curves[[i]][[2]],
            lty=lty[i], lwd=lwd[i], col=col[i], type='s')

  if(aehaz || length(times)) {
    un <- if(units == ' ' | units == '') '' else
      paste('/', tolower(units), sep='')
    haz <- round(nevents / totaltime, 4)
    txt <- paste(nevents, 'events')
    if(aehaz) txt <- paste(txt, ', hazard=', haz, un, sep='')
    if(length(times)) txt <- paste(txt, ', ', sep='')
    if(length(times)) txt <- paste(txt, cuminc)
    if(! labelc)
      text(xlim[2], ylim[2], txt, adj=1)
    else {
      maxlen <- max(nchar(sleva))
      sleva <- substring(paste(sleva, '                               '),
                         1, maxlen)
      for(j in 1 : ns)
        sleva[j] <- eval(parse(text=sprintf("expression(paste('%s   ',scriptstyle('(%s)')))", sleva[j], txt[j])))
    }
  }
  if(labelc) labcurve(curves, sleva, type='s', lty=lty, lwd=lwd,
                      opts=label.curves, col.=col)
  
  invisible(slev)
}


survdiffplot <-
  function(fit, order=1:2, fun=function(y) y, xlim, ylim, xlab,
           ylab="Difference in Survival Probability", time.inc,
           conf.int, conf=c("shaded", "bands", "diffbands", "none"),
           add=FALSE, lty=1, lwd=par('lwd'), col=1,
           n.risk=FALSE,  grid=NULL,
           srt.n.risk=0, adj.n.risk=1,
           y.n.risk, cex.n.risk=.6,
           cex.xlab=par('cex.lab'), cex.ylab=cex.xlab,
           convert=function(f) f) {

  conf <- match.arg(conf)
  if(missing(conf.int)) conf.int <- fit$conf.int
  if(! length(conf.int) | conf == "none") conf.int <- 0

  opar <- par(c('xpd', 'mar'))
  on.exit(par(opar))
  
  units <- Punits(fit$units, adds=FALSE, upfirst=TRUE, default='Day')
  maxtime <- fit$maxtime
  if(!length(maxtime)) maxtime <- max(fit$time)
  mintime <- min(fit$time, 0)
  pret    <- pretty(c(mintime, maxtime))
  maxtime <- max(pret)
  mintime <- min(pret)
  if(missing(time.inc)) {
    time.inc <- switch(units, Day=30, Month=1, Year=1, (maxtime - mintime) / 10)
    if(time.inc > maxtime) time.inc <- (maxtime - mintime) / 10
  }
  if(n.risk && !length(fit$n.risk)) {
    n.risk <- FALSE
    warning("fit does not have number at risk\nIs probably from a parametric model\nn.risk set to FALSE")
    }

  if(missing(xlab)) xlab <- if(units==' ') 'Time' else paste(units, "s", sep="")
  
  if(missing(xlim)) xlim <- c(mintime, maxtime)
  
  if(length(grid) && is.logical(grid)) grid <- if(grid) gray(.8) else NULL
  
  polyg <- ordGridFun(grid=FALSE)$polygon

  times <- sort(unique(c(fit$time, seq(mintime, maxtime, by=time.inc))))

    ## Note: summary.survfit computes standard errors on S(t) scale

  f <- convert(summary(fit, times=times))
  
  slev <- levels(f$strata)
  ns <- length(slev)
  if(ns !=2 ) stop('must have exactly two strata')
  a <- f$strata == slev[1]
  h <- function(level, times, f) {
    strata <- f$strata
    i     <- strata == level
    tim   <- f$time[i]
    surv  <- f$surv[i]
    se    <- f$std.err[i]
    nrisk <- f$n.risk[i]
    j     <- match(times, tim)
    list(time=times, surv=surv[j], se=se[j], nrisk=nrisk[j])
  }
  a <- h(slev[order[1]], times, f)
  b <- h(slev[order[2]], times, f)


  surv  <- if(conf == 'diffbands') (fun(a$surv) + fun(b$surv)) / 2
   else fun(a$surv) - fun(b$surv)
  se    <- sqrt(a$se^2 + b$se^2)

    z  <- qnorm((1 + conf.int) / 2)
  if(conf == 'diffbands') {
    lo <- surv - 0.5 * z * se
    hi <- surv + 0.5 * z * se
    k <- ! is.na(times + lo + hi) & times < xlim[2]
    polyg(c(times[k], rev(times[k])), c(lo[k], rev(hi[k])),
           col=gray(.9), type='s')
    return(invisible(slev))
  }
  lo <- surv - z * se
  hi <- surv + z * se

  if(missing(ylim)) ylim <- range(c(lo, hi), na.rm=TRUE)

  if(!add) {
    plot(times, surv, type='n', axes=FALSE,
         xlim=xlim, ylim=ylim, xlab='', ylab='')
    mgp.axis(2, labels=TRUE, axistitle=ylab, cex.lab=cex.ylab)
    mgp.axis(1, at=seq(xlim[1], max(pretty(xlim)), time.inc),
             labels=TRUE, axistitle=xlab, cex.lab=cex.xlab)
    }

  if(length(grid)) {
    xlm <- pretty(xlim)
    xlm <-   c(xlm[1], xlm[length(xlm)])
    xp  <- seq(xlm[1], xlm[2], by=time.inc)
    ylm <- pretty(ylim)
    yp  <- seq(min(ylm), max(ylm), by=ylm[2] - ylm[1])
    abline(h=yp, v=xp, col=grid, xpd=FALSE)
  }

  k <- !is.na(times + lo + hi)
  
  switch(conf,
         shaded=polyg(c(times[k], rev(times[k])), c(lo[k], rev(hi[k])),
           col=gray(.85), type='s'),
         bands={
           lines(times, lo, col=gray(.7))
           lines(times, hi, col=gray(.7))
         },
         diffbands=NULL,
         none=NULL)
  lines(times, surv, type='s', lwd=lwd, col=col)
  abline(h=0, col=gray(.7))
  title(sub=paste(slev[order], collapse=' - '), adj=0)

  if(n.risk) {
    nrisktimes <- seq(0, maxtime, by=time.inc)
    nriskinfo  <- convert(summary(fit, times=nrisktimes))
    anr        <- h(slev[order[1]], nrisktimes, nriskinfo)
    bnr        <- h(slev[order[2]], nrisktimes, nriskinfo)
    nrisk      <- pmin(anr$nrisk, bnr$nrisk)
    xd         <- xlim[2] - xlim[1]
    yd         <- ylim[2] - ylim[1]
    
    if(! add) {
      mar <- opar$mar
      if(mar[4] < 4) {mar[4] <- mar[4] + 2; par(mar=mar)}
    }
    par(xpd=NA)
    
    tt <- nrisktimes
    tt[1] <- xlim[1]
    if(missing(y.n.risk)) y.n.risk <- ylim[1]
    if(y.n.risk == 'auto') y.n.risk <- - diff(ylim) / 3
    yy  <- y.n.risk
    nri <- nrisk
    nri[tt > xlim[2]] <- NA
    text(tt[1], yy, nri[1], cex=cex.n.risk,
         adj=adj.n.risk, srt=srt.n.risk)
    if (length(nri) > 1)
      text(tt[-1], yy, nri[-1], cex=cex.n.risk, adj=1)
  }
  invisible(slev)
}
