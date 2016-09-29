survplot <- function(fit, ...) UseMethod("survplot")

survplot.rms <-
  function(fit, ..., xlim, 
           ylim=if(loglog) c(-5,1.5) else 
            if(what == "survival" & missing(fun)) c(0,1),
           xlab, ylab, time.inc,
           what=c("survival","hazard"),
           type=c("tsiatis","kaplan-meier"),
           conf.type=c("log","log-log","plain","none"),
           conf.int=FALSE, conf=c("bands","bars"), mylim=NULL,
           add=FALSE, label.curves=TRUE,
           abbrev.label=FALSE, levels.only=FALSE,
           lty, lwd=par('lwd'),
           col=1, col.fill=gray(seq(.95, .75, length=5)),
           adj.subtitle=TRUE, loglog=FALSE, fun, n.risk=FALSE, logt=FALSE,
           dots=FALSE, dotsize=.003, grid=NULL,
           srt.n.risk=0, sep.n.risk=.056, adj.n.risk=1,
           y.n.risk, cex.n.risk=.6, cex.xlab=par('cex.lab'),
           cex.ylab=cex.xlab, pr=FALSE)
{

  what <- match.arg(what)
  polyg <- ordGridFun(grid=FALSE)$polygon
  ylim <- ylim  ## before R changes missing(fun)
  type <- match.arg(type)
  conf.type <- match.arg(conf.type)
  conf <- match.arg(conf)

  opar <- par(c('mar', 'xpd'))
  on.exit(par(opar))

  cylim <- function(ylim)
    if(length(mylim)) c(min(ylim[1], mylim[1]), max(ylim[2], mylim[2]))
    else ylim

  psmfit <- inherits(fit,'psm')
  if(what == "hazard" && !psmfit)
    stop('what="hazard" may only be used for fits from psm')
  if(what == "hazard" & conf.int > 0) {
    warning('conf.int may only be used with what="survival"')
    conf.int <- FALSE
  }
  
  if(loglog) {
    fun <- function(x) logb(-logb(ifelse(x == 0 | x == 1, NA, x)))
    use.fun <- TRUE
  }
  else if(!missing(fun)) {
    use.fun <- TRUE
    if(loglog) stop("cannot specify loglog=T with fun")
  }
  else {
    fun <- function(x) x
    use.fun <- FALSE
  }
  
  if(what == "hazard" & loglog)
    stop('may not specify loglog=T with what="hazard"')

  if(use.fun | logt | what == "hazard") { dots <- FALSE; grid <- NULL }
  
  cox <- inherits(fit,"cph")
  if(cox) {
    if(n.risk | conf.int>0) surv.sum <- fit$surv.summary
    exactci <- !(is.null(fit$x)|is.null(fit$y))
    ltype <- "s"	#step functions for cph
  }
  
  else {
    if(n.risk) stop("the n.risk option applies only to fits from cph")
    exactci <- TRUE
    ltype <- "l"
  }
  
  par(xpd=NA)
  
  ## Compute confidence limits for survival based on -log survival,
  ## constraining to be in [0,1]; d = std.error of cum hazard * z value
  ciupper <- function(surv, d) ifelse(surv == 0, 0, pmin(1, surv*exp(d)))
  cilower <- function(surv, d) ifelse(surv == 0, 0, surv*exp(-d))
  
  labelc <- is.list(label.curves) || label.curves
  
  units <- fit$units
  if(missing(ylab)) {
    if(loglog) ylab <- "log(-log Survival Probability)"
    else if(use.fun) ylab <- ""
    else if(what == "hazard") ylab <- "Hazard Function"
    else ylab <- "Survival Probability"
  }
  if(missing(xlab)) {
    if(logt) xlab <- paste("log Survival Time in ",units,"s",sep="")
    else xlab <- paste(units,"s",sep="")
  }
  
  maxtime <- fit$maxtime
  maxtime <- max(pretty(c(0,maxtime)))
  if(missing(time.inc)) time.inc <- fit$time.inc

  if(missing(xlim)) 
    xlim <- if(logt)logb(c(maxtime/100,maxtime)) else c(0,maxtime)

  if(length(grid) && is.logical(grid)) grid <- if(grid) gray(.8) else NULL

  if(is.logical(conf.int)) {
    if(conf.int) conf.int <- .95	else conf.int <- 0
  }
  zcrit <- qnorm((1+conf.int)/2)
  
  xadj <- Predict(fit, type='model.frame', np=5,
                  factors=rmsArgs(substitute(list(...))))
  info <- attr(xadj, 'info')
  varying <- info$varying
  if(length(varying) > 1)
    stop('cannot vary more than one predictor')

  adjust <- if(adj.subtitle) info$adjust else NULL
  if(length(xadj)) {
    nc <- nrow(xadj)
    covpres <- TRUE
  }
  else {
    nc <- 1
    covpres <- FALSE
  }
  y <- if(length(varying)) xadj[[varying]] else ''

  curve.labels <- NULL
  xd <- xlim[2] - xlim[1]
  if(n.risk & !add) {
    mar <- opar$mar
    if(mar[4] < 4) {
      mar[4] <- mar[4] + 2
      par(mar=mar)
    }
  }
  
  ## One curve for each value of y, excl style used for C.L.
  lty <- if(missing(lty)) seq(nc+1)[-2]
  else
    rep(lty, length=nc)
  col <- rep(col, length=nc)
  lwd <- rep(lwd, length=nc)
  
  i <- 0
  if(levels.only) y <- gsub('.*=', '', y)
  abbrevy <- if(abbrev.label) abbreviate(y) else y
  abbrevy <- if(is.factor(abbrevy)) as.character(abbrevy) else format(abbrevy)
  
  if(labelc || conf == 'bands') curves <- vector('list',nc)

  for(i in 1:nc) {
    ci <- conf.int
    ay <- if(length(varying)) xadj[[varying]] else ''
    if(covpres) {
      adj <- xadj[i,,drop=FALSE]
      w <- survest(fit, newdata=adj, fun=fun, what=what,
                   conf.int=ci, type=type, conf.type=conf.type)
    }
    else w <- survest(fit, fun=fun, what=what, conf.int=ci,
                      type=type, conf.type=conf.type)

    time <- w$time
    if(logt) time <- logb(time)
    s <- !is.na(time) & (time>=xlim[1])
    surv <- w$surv
    if(is.null(ylim)) ylim <- cylim(range(surv, na.rm=TRUE))
    stratum <- w$strata
    if(is.null(stratum)) stratum <- 1
    if(!is.na(stratum)) {
      ##can be NA if illegal strata combinations requested
      cl <- if(is.factor(ay)) as.character(ay)
      else format(ay)
      curve.labels <- c(curve.labels, abbrevy[i])
      if(i == 1 & !add) {				
        plot(time, surv, xlab='', xlim=xlim,
             ylab='', ylim=ylim, type="n", axes=FALSE)	
        mgp.axis(1, at=if(logt)pretty(xlim) else
                 seq(xlim[1], max(pretty(xlim)), time.inc), labels=TRUE,
                 axistitle=xlab, cex.lab=cex.xlab)
        mgp.axis(2, at=pretty(ylim), axistitle=ylab, cex.lab=cex.ylab)
        
        if(!logt & (dots || length(grid))) {
          xlm <- pretty(xlim)
          xlm <-  c(xlm[1], xlm[length(xlm)])
          xp <- seq(xlm[1], xlm[2],by=time.inc)
          yd <- ylim[2] - ylim[1]
          if(yd <= .1) yi <- .01
          else if(yd <= .2) yi <- .025
          else if(yd <= .4) yi <- .05
          else yi <- .1
          yp <- seq(ylim[2], ylim[1] +
                    if(n.risk && missing(y.n.risk)) yi else 0, 
                    by=- yi)
          if(dots) for(tt in xp)symbols(rep(tt, length(yp)), yp,
                                        circles=rep(dotsize, length(yp)),
                                        inches=dotsize, add=TRUE)
          else abline(h=yp, v=xp, col=grid, xpd=FALSE)
        }
      }
      tim <- time[s]; srv <- surv[s]
      if(conf.int > 0 && conf == 'bands') {
        blower <- w$lower[s]
        bupper <- w$upper[s]
      }
      if(max(tim) > xlim[2]) {
        if(ltype == "s") {
          ##Get estimate at last permissible point to plot
          ## s.last <- min(srv[tim<=xlim[2]+1e-6])  #not work with function
          s.last <- srv[tim <= xlim[2] + 1e-6]
          s.last <- s.last[length(s.last)]
          k <- tim < xlim[2]
          tim <- c(tim[k], xlim[2]); srv <- c(srv[k], s.last)
          if(conf.int > 0 && conf == 'bands') {
            low.last <- blower[time <= xlim[2] + 1e-6]
            low.last <- low.last[length(low.last)]
            up.last  <- bupper[time <= xlim[2] + 1e-6]
            up.last  <- up.last[length(up.last)]
            blower   <- c(blower[k],low.last)
            bupper   <- c(bupper[k],up.last)
          }
        }
        else tim[tim > xlim[2]] <- NA
      }
      
      ##don't let step function go beyond x-axis -
      ##this cuts it off but allows step to proceed axis end
      
      if(conf != 'bands')
        lines(tim, srv, type=ltype, lty=lty[i], col=col[i], lwd=lwd[i])
      
      if(labelc || conf == 'bands') curves[[i]] <- list(tim, srv)
      
      if(pr) {
        zest <- rbind(tim,srv)
        dimnames(zest) <- list(c("Time","Survival"), rep("",length(srv)))
        cat("\nEstimates for ", cl,"\n\n")
        print(zest, digits=3)
      }
      if(conf.int > 0) {
        if(conf == "bands") {
          polyg(x = c(tim,    rev(tim)),
                y = c(blower, rev(bupper)),
                col =  col.fill[i], type=ltype)
        }
        else {
          if(exactci) { # not from cph(surv=T)
            tt <- seq(0,  maxtime, time.inc)
            v <- survest(fit, newdata=adj, times=tt,
                         what=what, fun=fun,
                         conf.int=ci, type=type, conf.type=conf.type)
            tt    <- v$time   #may not get predictions at all t
            ss    <- v$surv
            lower <- v$lower
            upper <- v$upper
            if(!length(ylim)) ylim <- cylim(range(ss, na.rm=TRUE))
            if(logt) tt <- logb(ifelse(tt == 0, NA, tt))
          }
          else {
            tt <- as.numeric(dimnames(surv.sum)[[1]])
            if(logt) tt <- logb(tt)
            ss <- surv.sum[,stratum,'Survival']^
            exp(w$linear.predictors)
            se <- surv.sum[,stratum,'std.err']
            ss <- fun(ss)
            lower <- fun(cilower(ss, zcrit*se))
            upper <- fun(ciupper(ss, zcrit*se))
            ss[is.infinite(ss)] <- NA
            lower[is.infinite(lower)] <- NA
            upper[is.infinite(upper)] <- NA
          }
          tt <- tt + xd*(i-1)*.01
          errbar(tt, ss, upper, lower, add=TRUE, lty=lty[i], col=col[i])
        }
      }
      if(n.risk) {
        if(length(Y <- fit$y)) {
          tt <- seq(max(0,xlim[1]),min(maxtime,xlim[2]),by=time.inc)
          ny <- ncol(Y)
          if(!length(str <- fit$strata)) Y <- Y[,ny-1]
          else Y <- Y[unclass(str) == unclass(stratum), ny - 1]
          nrisk <-
            rev(cumsum(table(
              cut(-Y,sort(unique(-c(tt,range(Y)+
                                    c(-1,1)))))
              )[-length(tt)-1]))	
        }
        else {
          if(is.null(surv.sum))
            stop("you must use surv=T or y=T in fit to use n.risk=T")
          tt <- as.numeric(dimnames(surv.sum)[[1]])
          l <- (tt >= xlim[1]) & (tt <= xlim[2])
          tt <- tt[l]
          nrisk <- surv.sum[l,stratum,2]
        }
        tt[1] <- xlim[1]  #was xd*.015, .030, .035
        yd <- ylim[2] - ylim[1]
        if(missing(y.n.risk)) y.n.risk <- ylim[1]
        if(y.n.risk == 'auto') y.n.risk <- - diff(ylim) / 3
        yy <- y.n.risk + yd*(nc-i)*sep.n.risk #was .029, .038, .049
        
        nri <- nrisk
        nri[tt > xlim[2]] <- NA
        text(tt[1], yy, nri[1], cex=cex.n.risk,
             adj=adj.n.risk, srt=srt.n.risk)
        if (length(nri) > 1)
          text(tt[-1], yy, nri[-1], cex=cex.n.risk, adj=1)
        text(xlim[2]+xd*.025, yy, adj=0, curve.labels[i], cex=cex.n.risk)
      }
    }
  }
  
  ## to keep bands from covering up lines plot lines last
  if(conf == 'bands') for(i in 1:length(y))
    lines(curves[[i]][[1]], curves[[i]][[2]], type=ltype, lty=lty[i],
          col=col[i], lwd=lwd[i])

  if(labelc) labcurve(curves, curve.labels, type=ltype, lty=lty, col.=col,
                      lwd=lwd, opts=label.curves)

  if(length(adjust)) title(sub=paste("Adjusted to:",adjust),
                           adj=0, cex=.6)
  
  invisible(list(adjust=adjust, curve.labels=curve.labels))
}
