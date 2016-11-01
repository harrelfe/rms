#Compute various measures of reliability and discrimination for a set
#of predicted probabilities p or predicted logits logit.
#If pl=T, the following apply:
#  Plots reliability curve, for which xlab is optional label.
#  If smooth=T and pl=T, plots lowess(p,y,iter=0)
#  lim is x-axis and y-axis range, default=c(0,1)
#  If m or g is specified, also computes and plots proportions of y=1
#  by quantile groups of p (or 1/(1+exp(-logit))).  If m is given,
#  groups are constructed to have m observations each on the average.
#  Otherwise, if g is given, g quantile groups will be constructed.
#  If instead cuts is given, proportions will be computed based on the
#  cut points in the vector cuts, e.g. cuts<-seq(0,1,by=.2). 
#  If legendloc is given, a legend will be plotted there
#  Otherwise, it is placed at (.6, .38)
#  Use legendloc=locator(1) to use the mouse for legend positioning.
#  Use legendloc="none" to suppress legend.
#  If statloc is given, some statistics will be plotted there
#  Use statloc=locator(1) to use the mouse.  This is done after the legend.
#  legendloc and statloc can be lists as returned by locator() or they
#  can be vectors, e.g. c(x,y).
#
#Frank Harrell 1 Jun 91
#
val.prob <- function(p, y, logit, group, weights=rep(1,length(y)),
                     normwt=FALSE, pl=TRUE, smooth=TRUE, logistic.cal=TRUE,
					 xlab="Predicted Probability", ylab="Actual Probability",
					 lim=c(0,1), m, g, cuts, emax.lim=c(0,1),
					 legendloc=lim[1] + c(.55 * diff(lim), .27 * diff(lim)),
					 statloc=c(0,.99), riskdist=c("predicted", "calibrated"),
           cex=.7, mkh=.02,
					 connect.group=FALSE, connect.smooth=TRUE, 
					 g.group=4, evaluate=100, nmin=0)
{

  if(missing(p)) p <- plogis(logit)	else logit <- qlogis(p)
  if(length(p) != length(y)) stop("lengths of p or logit and y do not agree")
  names(p) <- names(y) <- names(logit) <- NULL

  riskdist <- match.arg(riskdist)

  Spi <- function(p, y) {
    z <- sum((y - p)*(1 - 2*p)) /
      sqrt(sum((1 - 2 * p) * (1 - 2 * p) * p * (1-p)))
    P <- 2 * pnorm(- abs(z))
    c(Z=z, P=P)
  }
  
  if(! missing(group)) {
    if(length(group)==1 && is.logical(group) && group)
      group <- rep('', length(y))
    if(! is.factor(group)) group <- 
      if(is.logical(group) || is.character(group)) 
        as.factor(group) else cut2(group, g=g.group)
    names(group) <- NULL
    nma <- ! (is.na(p + y + weights) | is.na(group))
    ng  <- length(levels(group))
  }
  else {
    nma <- ! is.na(p + y + weights)
    ng <- 0
  }
  
  logit <- logit[nma]
  y     <- y[nma]
  p     <- p[nma]
  if(ng > 0) {
    group   <- group[nma]
    weights <- weights[nma]
    return(val.probg(p, y, group, evaluate, weights, normwt, nmin) )
  }
  
  if(length(unique(p)) == 1) {
    P     <- mean(y)
    Intc  <- qlogis(P)
    n     <- length(y)
    D     <- -1 / n
    L01   <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
    L.cal <- -2 * sum(y * Intc  - logb(1 + exp(Intc)),  na.rm=TRUE)
    U.chisq <- L01 - L.cal
    U.p <- 1 - pchisq(U.chisq, 1)
    U <- (U.chisq - 1) / n
    Q <- D - U
    spi <- unname(Spi(p, y))
    stats <- c(0, .5, 0, D, 0, 1, U, U.chisq, U.p, Q,
               mean((y - p[1]) ^ 2), Intc, 0, 0, 0,
               rep(abs(p[1] - P), 2), spi)
    names(stats) <- c("Dxy","C (ROC)", 
                      "R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
                      "Brier","Intercept","Slope","Emax","E90","Eavg",
                      "S:z", "S:p")
    return(stats)
  }
  
  i <- ! is.infinite(logit)
  nm <- sum(! i)
  if(nm > 0) warning(paste(nm,
    "observations deleted from logistic calibration due to probs. of 0 or 1"))
  f.fixed <- lrm.fit(logit[i], y[i], initial=c(0., 1.), maxit=1L)
  f.recal <- lrm.fit(logit[i], y[i])
  stats <- f.fixed$stats
  n <- stats["Obs"]
  predprob <- seq(emax.lim[1], emax.lim[2], by=.0005)

  Sm <- lowess(p, y, iter=0)
  cal.smooth <- approx(Sm, xout=p, ties=mean)$y
  er   <- abs(p - cal.smooth)
  eavg <- mean(er)
  emax <- max(er)
  e90  <- unname(quantile(er, 0.9))
  
  if(pl) {
    plot(.5, .5, xlim=lim, ylim=lim, type="n", xlab=xlab, ylab=ylab)
    abline(0, 1, lwd=6, col=gray(.85))
    lt <- 1; leg <- "Ideal"; marks <- -1; lwd <- 6; col <- gray(.85)
    if(logistic.cal) {
      lt <- c(lt, 1); leg <- c(leg, "Logistic calibration")
      lwd <- c(lwd, 1); col <- c(col, 'black')
      marks <- c(marks, -1)
    }
    if(smooth) {
      if(connect.smooth) { 
        lines(Sm, lty=3)
        lt <- c(lt, 3)
        lwd <- c(lwd, 1); col <- c(col, 'black')
        marks <- c(marks, -1)
      }
      else {
        points(Sm)
        lt <- c(lt, 0)
        lwd <- c(lwd, 1); col <- c(col, 'black')
        marks <- c(marks, 1)
      }
      leg <- c(leg, "Nonparametric")
    }
    if(! missing(m) | ! missing(g) | ! missing(cuts)) {
           if(! missing(m))    q <- cut2(p, m=m, levels.mean=TRUE, digits=7)
      else if(! missing(g))    q <- cut2(p, g=g, levels.mean=TRUE, digits=7)
      else if(! missing(cuts)) q <- cut2(p, cuts=cuts, levels.mean=TRUE,
                                        digits=7)
      means <- as.numeric(levels(q))
      prop <- tapply(y, q, function(x) mean(x, na.rm=TRUE))
      points(means, prop, pch=2)
      if(connect.group) {lines(means, prop); lt <- c(lt, 1)}
      else lt <- c(lt, 0)
           leg <- c(leg, "Grouped observations")
      col <- c(col, 'black'); lwd <- c(lwd, 1)     
      marks <- c(marks, 2)
    }
	}
  lr <- stats["Model L.R."]
  p.lr <- stats["P"]
  D <- (lr - 1) / n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
  U.chisq <- L01 - f.recal$deviance[2]
  p.U <- 1 - pchisq(U.chisq, 2)
  U <- (U.chisq - 2)/n
  Q <- D - U
  Dxy <- stats["Dxy"]
  C   <- stats["C"]
  R2  <- stats["R2"]
  B   <- mean((p - y) ^ 2)
  spi <- unname(Spi(p, y))
  stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B, f.recal$coef,
             emax, e90, eavg, spi)
  names(stats) <- c("Dxy","C (ROC)", 
                    "R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
                    "Brier","Intercept","Slope","Emax","E90","Eavg","S:z","S:p")

  if(pl) {
    logit <- seq(-7, 7, length=200)
    prob  <- plogis(logit)
    pred.prob <- f.recal$coef[1] + f.recal$coef[2] * logit
    pred.prob <- plogis(pred.prob)
    if(logistic.cal) lines(prob, pred.prob, lty=1)
    lp <- legendloc
    if(! is.logical(lp)) {
      if(! is.list(lp)) lp <- list(x=lp[1],y=lp[2])
      legend(lp, leg, lty=lt, pch=marks, cex=cex, lwd=lwd, col=col, bty="n")
    }
    if(! is.logical(statloc)) {
      dostats <- c("Dxy", "C (ROC)", "R2", "D", "U", "Q", "Brier",
                   "Intercept", "Slope", "Emax", "E90", "Eavg", 
                   "S:z", "S:p")      
  
      leg <- format(names(stats)[dostats]) #constant length
      leg <- paste(leg, ":", format(stats[dostats]),sep="")
      if(! is.list(statloc)) statloc <- list(x=statloc[1], y=statloc[2])
      text(statloc,paste(format(names(stats[dostats])), collapse="\n"),
           adj=c(0, 1), cex=cex)
      text(statloc$x + .225 * diff(lim), statloc$y,
           paste(format(round(stats[dostats], 3)),
                 collapse="\n"), adj=c(1,1), cex=cex)
    }
    if(is.character(riskdist)) {
      if(riskdist=="calibrated") {
        x <- f.recal$coef[1] + f.recal$coef[2] * qlogis(p)
        x <- plogis(x)
        x[p == 0] <- 0; x[p == 1] <- 1
      } else x <- p
      bins <- seq(lim[1], lim[2], length=101)
      x <- x[x >= lim[1] & x <= lim[2]]
      f <- table(cut(x, bins))
      j <- f > 0
      bins <- (bins[-101])[j]
      f <- f[j]
      f <- lim[1] + .15 * diff(lim) * f / max(f)
      segments(bins, 0, bins, f)
    }
  }	
  stats
}


val.probg <- function(p, y, group, evaluate=100, weights, normwt, nmin)
{
  if(normwt) weights <- length(y)*weights/sum(weights)
  ng <- length(lg <- levels(group))
  if(ng==1) {ng <- 0; lg <- character(0)}
  stats <- matrix(NA, nrow=ng+1, ncol=12,
				  dimnames=list(nn <- c(lg,'Overall'), 
					c('n','Pavg','Obs','ChiSq','ChiSq2','Eavg',
					  'Eavg/P90','Med OR','C','B','B ChiSq','B cal')))
  curves <- vector('list',ng+1)
  names(curves) <- nn
  q.limits <- c(.01,.025,.05,.1,.25,.5,.75,.9,.95,.975,.99)
  limits <- matrix(NA, nrow=ng+1, ncol=length(q.limits),
				   dimnames=list(nn, as.character(q.limits)))
  for(i in 1:(ng+1))
    {
      s <- if(i==(ng+1)) 1:length(p) else group==lg[i]
      P <- p[s]
      Y <- y[s]
      wt <- weights[s]
      lims <- wtd.quantile(P, wt, q.limits, na.rm=FALSE, normwt=FALSE)
      limits[i,] <- lims
      n <- sum(wt)
      n1 <- sum(wt[Y == 1])
      c.index <- (mean(wtd.rank(P, wt, na.rm=FALSE, normwt=FALSE)[Y == 1]) - 
                  (n1 + 1)/2)/(n - n1)
      ## c.index <- somers2(P, Y, wt, normwt=FALSE, na.rm=FALSE)['C']
      sm <- wtd.loess.noiter(P, Y, wt, na.rm=FALSE, type='all')  
      ##all -> return all points
      curve <- if(length(sm$x) > evaluate)
        approx(sm, xout=seq(min(P), max(P), length=evaluate), ties=mean) else
      {
        o <- order(sm$x)
        nd <- ! duplicated(sm$x[o])
        list(x=sm$x[o][nd], y=sm$y[o][nd])
      }
      if(nmin > 0)
        {
          cuts <- wtd.quantile(P, wt, c(nmin, n-nmin)/n, normwt=FALSE, na.rm=FALSE)
          keep <- curve$x >= cuts[1] & curve$x <= cuts[2]
          curve <- list(x=curve$x[keep], y=curve$y[keep])
        }
      curves[[i]] <- curve
      cal.smooth <- sm$y
      eavg <- sum(wt * abs(P - cal.smooth))/n
      b    <- sum(wt * ((P - Y)^2))/n
      E0b  <- sum(wt * P * (1 - P))/n
      Vb   <- sum(wt * ((1 - 2 * P)^2) * P * (1 - P))/n/n
      bchisq <- (b - E0b)^2 / Vb
      b.cal  <- sum(wt * ((cal.smooth - Y)^2))/n
      
      pred  <- sum(wt * P)/n
      obs   <- sum(wt * Y)/n
      L <- ifelse(P==0 | P==1, NA, qlogis(P))
      w <- ! is.na(L)
      del <- matrix(c(sum((wt*(Y-P))[w]),sum((wt*L*(Y-P))[w])),ncol=2)
      v <- rbind(c(sum((wt*P*(1-P))[w]), sum((wt*L*P*(1-P))[w])),
                 c(NA, sum((wt*L*L*P*(1-P))[w])))
      v[2,1] <- v[1,2]
      chisq  <- (sum(wt * (P - Y))^2) / sum(wt * P * (1 - P))
      chisq2 <- del %*% solve(v) %*% t(del)
      p90    <- diff(lims[c(3,9)])
      Lcal   <- ifelse(cal.smooth <= 0 | cal.smooth >= 1, NA,
                       qlogis(cal.smooth))
      or <- exp(wtd.quantile(abs(L - Lcal), wt, .5, na.rm=TRUE, normwt=FALSE))
      stats[i,] <- c(n, pred, obs, chisq, chisq2, eavg, eavg/p90, or, c.index,
                     b, bchisq, b.cal)
    }
  structure(list(stats=stats, cal.curves=curves, quantiles=limits), 
			class='val.prob')
}

print.val.prob <- function(x, ...)
{
  print(round(x$stats,3))
  cat('\nQuantiles of Predicted Probabilities\n\n')
  print(round(x$quantiles,3))
  invisible()
}

plot.val.prob <- function(x, 
						  xlab="Predicted Probability", 
						  ylab="Actual Probability",
						  lim=c(0,1), statloc=lim, stats=1:12, cex=.5, 
						  lwd.overall=4, quantiles=c(0.05,0.95),
						  flag=function(stats) ifelse(
						   stats[,'ChiSq2'] > qchisq(.99,2) |
						   stats[,'B ChiSq'] > qchisq(.99,1),'*',' '), ...)
{
  stats <- x$stats[,stats,drop=FALSE]
  lwd <- rep(par('lwd'), nrow(stats))
  lwd[dimnames(stats)[[1]]=='Overall'] <- lwd.overall
  curves <- x$cal.curves

  labcurve(curves, pl=TRUE, xlim=lim, ylim=lim, 
		   xlab=xlab, ylab=ylab, cex=cex, lwd=lwd, ...)
  abline(a=0, b=1, lwd=6, col=gray(.86))
  if(is.logical(statloc) && ! statloc) return(invisible())

  if(length(quantiles))
    {
      limits <- x$quantiles
      quant <- round(as.numeric(dimnames(limits)[[2]]),3)
      w <- quant %in% round(quantiles,3)
      if(any(w)) for(j in 1:nrow(limits))
        {
          qu <- limits[j,w]
          scat1d(qu, y=approx(curves[[j]], xout=qu, ties=mean)$y)
        }
    }

  xx <- statloc[1]; y <- statloc[2]
  for(i in 0:ncol(stats))
    {
      column.text <- if(i==0) c('Group',
                          paste(flag(stats),dimnames(stats)[[1]],sep='')) else
	  c(dimnames(stats)[[2]][i], 
		format(round(stats[,i], if(i %in% c(4:5,11))1 else 3)))
      cat(column.text, '\n')
      text(xx, y, paste(column.text, collapse='\n'), adj=0, cex=cex)
      xx <- xx + (1 + .8 * max(nchar(column.text))) * cex * par('cxy')[1]
    }
  invisible()
}
