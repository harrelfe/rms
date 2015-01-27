residuals.lrm <-
  function(object, 
           type=c("li.shepherd", "ordinary","score","score.binary","pearson",
             "deviance","pseudo.dep","partial",
             "dfbeta","dfbetas","dffit","dffits","hat","gof","lp1"),
           pl=FALSE, xlim, ylim, kint=1, label.curves=TRUE, 
           which, ...)
{
  gotsupsmu <- FALSE
  type <- match.arg(type)
  dopl <- (is.logical(pl) && pl) || is.character(pl)
  ylabpr <- NULL   # y-axis label for partial residuals

  k <- object$non.slopes
  L <- object$linear.predictors
  isorm <- inherits(object, 'orm')

  trans <- object$trans
  cumprob <- if(length(trans)) trans$cumprob else plogis
  if(length(L)==0) stop('you did not use linear.predictors=TRUE for the fit')

  if(kint < 1 | kint > k) stop(paste('kint must be from 1-', k, sep=''))

  cof <- object$coefficients
  ordone <- type %in% c('li.shepherd','partial','gof','score','score.binary')
  ## residuals explicitly handled for ordinal model
  if(ordone && !missing(kint)) 
    stop('may not specify kint for li.shepherd, partial, score, score.binary, or gof')

  if(isorm) L <- L - cof[attr(L, 'intercepts')] + cof[1]

  if(k > 1 && kint !=1 && !ordone) L <- L - cof[1] + cof[kint]
  P <- cumprob(L)
  if(length(Y <- object$y) == 0) stop("you did not specify y=TRUE in the fit")
  rnam <- names(Y)
  cnam <- names(cof)
  if(!is.factor(Y)) Y <- factor(Y)
  lev  <- levels(Y)
  lev2 <- names(cof)[1:k]
  Y <- unclass(Y) - 1
  if(!ordone && k > 1) Y <- Y >= kint
  if(k > 1 && missing(kint) && !ordone)
    warning(paste('using first intercept and ',
                  lev2[kint],
                  ' to compute residuals or test GOF',
                  sep=''))
  
  if(type=="gof") {
    if(length(X <- object$x)==0)
      stop("you did not use x=TRUE in the fit")
    stats <- matrix(NA, nrow=k, ncol=5,
                    dimnames=list(if(k > 1) lev2,
                      c("Sum of squared errors", "Expected value|H0",
                        "SD", "Z", "P")))
    X <- cbind(1,X)
    for(j in 1:k) {
      y <- Y >= j
      p <- cumprob(L - cof[1] + cof[j])
      sse <- sum((y - p)^2)
      wt <- p * (1 - p)
      d <- 1 - 2 * p
      z <- lm.wfit(X, d, wt, method='qr')
      ##    res <- summary(lm.wfit(X, d, wt, method="qr"))$residuals  11Apr02
      res <- z$residuals * sqrt(z$weights)
      sd <- sqrt(sum(res^2))
      ev <- sum(wt)
      z <- (sse-ev)/sd
      P <- 2 * (1 - pnorm(abs(z)))
      stats[j,] <- c(sse, ev, sd, z, P)
    }
    return(drop(stats))
  }
  
  naa <- object$na.action
  
  if(type=="ordinary") return(naresid(naa, Y - cumprob(L)))

  if(type %in% c('score','score.binary','partial')) {
    nc <- length(cof)
    if(missing(which)) which <- if(type == 'score') 1:nc else 1:(nc - k) else
    if(type=='score') which <- k + which
  }
  
  if(type=='score' || type=='score.binary')
    plotit <- function(w, ylim, xlab, ylab, lev=names(w)) {
      statsum <- function(x) {
        n <- length(x)
        xbar <- sum(x) / n
        if(n < 2) {low <- hi <- NA} else {
          se   <- 1.959964 * sqrt(sum((x - xbar)^2) / (n - 1) / n)
          low  <- xbar - se; hi <- xbar + se
        }
        c(mean=xbar, lower=low, upper=hi)
      }
      k <- length(w)
      w <- lapply(w, statsum)
      plot(c(1,k), c(0,0), xlab=xlab, ylab=ylab,
           ylim=if(length(ylim)==0) range(unlist(w)) else ylim,
           type='n', axes=FALSE)
      mgp.axis(2)
      mgp.axis(1, at=1:k, labels=lev)
      abline(h=0, lty=2, lwd=1)
      ii <- 0
      for(ww in w) {
        ii <- ii+1
        points(ii, ww[1])
        errbar(ii, ww[1], ww[3], ww[2], add=TRUE)
      }
    }
  
  if(type=='score.binary') {
    if(k==1)  stop('score.binary only applies to ordinal models')
    if(!dopl) stop('score.binary only applies if you are plotting')
    if(!length(X <- unclass(object$x)))
      stop('you did not specify x=TRUE for the fit')
    xname <- dimnames(X)[[2]]
    yname <- as.character(formula(object))[2]
    for(i in which) {
      xi <- X[,i]
      r <- vector('list',k)
      names(r) <- lev[-1]
      for(j in 1:k) 
        r[[j]] <- xi * ((Y >= j) - cumprob(L - cof[1] + cof[j]))
      if(pl!='boxplot') plotit(r, ylim=if(missing(ylim)) NULL else ylim,
           xlab=yname, ylab=xname[i])
      else
        boxplot(r, varwidth=TRUE, notch=TRUE, err=-1,
                ylim=if(missing(ylim)) quantile(unlist(r), c(.1, .9))
                else ylim, ...)
      title(xlab=yname, ylab=xname[i])
    }
    invisible()
  }
  
  if(type=="score") {
    if(!length(X <- unclass(object$x)))
      stop("you did not specify x=TRUE for the fit")
    if(isorm && object$family != 'logistic')
      stop('score residuals not yet implemented for orm with non-logistic family')
    if(k == 1) return(naresid(naa, cbind(1, X) * (Y - P))) # only one intercept
    z <- function(i, k, L, coef) cumprob(coef[pmin(pmax(i, 1), k)] + L)
    ## Mainly need the pmax - 0 subscript will drop element from vector
    ##  z$k <- k; z$L <- L-cof[1]; z$coef <- cof
    formals(z) <- list(i=NA, k=k, L=L-cof[1], coef=cof)
    ## set defaults in fctn def'n
    u <- matrix(NA, nrow=length(L), ncol=length(which),
                dimnames=list(names(L), names(cof)[which]))
    pq <- function(x) x * (1 - x)
    ## Compute probabilities of being in observed cells
    pc <- ifelse(Y==0, 1 - z(1), ifelse(Y == k, z(k), z(Y) - z(Y + 1)) )
    xname <- dimnames(X)[[2]]
    yname <- as.character(formula(object))[2]
    ii <- 0
    for(i in which) {
      ii <- ii + 1
      di  <- if(i <= k) ifelse(Y==0, if(i==1) 1 else 0, Y==i) else X[,i - k]
      di1 <- if(i <= k) ifelse(Y==0 | Y==k, 0, (Y + 1) == i)  else X[,i - k]
      ui  <- ifelse(Y == 0, -z(1) * di,
                    ifelse(Y == k, (1 - z(k)) * di,
                           (pq(z(Y)) * di - pq(z(Y + 1)) * di1) / pc ) )
      u[,ii] <- ui
      if(dopl && i > k) {
        if(pl=='boxplot') {
          boxplot(split(ui, Y), varwidth=TRUE, notch=TRUE,
                  names=lev, err=-1, 
                  ylim=if(missing(ylim)) quantile(ui, c(.1, .9))
                  else ylim, ...)
          title(xlab=yname, ylab=paste('Score Residual for', xname[i - k]))
        }
        else plotit(split(ui,Y), ylim=if(missing(ylim)) NULL else ylim,
                    lev=lev,
                    xlab=yname,
                    ylab=paste('Score Residual for', xname[i-k]))
      }
    }
    return(if(dopl) invisible(naresid(naa, u)) else naresid(naa, u))
  }
  if(type == "li.shepherd") {
    if(length(X <- object$x)==0)
      stop("you did not use x=TRUE in the fit")
    N <- length(Y)
    px <- 1 - cumprob(outer(cof[1:k],
                           as.vector(X %*% cof[- (1:k)]), "+"))
    low.x = rbind(0, px)[cbind(Y + 1L, 1:N)]
    hi.x  = 1 - rbind(px, 1)[cbind(Y + 1L, 1:N)]
    return(low.x - hi.x)
  }
  
  if(type=="pearson") return(naresid(naa, (Y - P) / sqrt(P * (1 - P))))
  
  if(type=="deviance") {
    r <- ifelse(Y==0,-sqrt(2 * abs(log(1 - P))), sqrt(2 * abs(logb(P))))
    return(naresid(naa, r))
  }
  
  if(type=="pseudo.dep") {
    r <- L + (Y - P) / P / (1 - P)
    return(naresid(naa, r))
  }
  
  if(type=="partial") {
    if(!length(X <- unclass(object$x)))
      stop("you did not specify x=TRUE in the fit")
    cof.int <- cof[1 : k]
    cof     <- cof[- (1 : k)]
    if(!missing(which)) {
      X <- X[, which, drop=FALSE]
      cof <- cof[which]
    }
    atx <- attributes(X)
    dx <- atx$dim
    if(k==1) r <- X * matrix(cof, nrow=dx[1], ncol=dx[2], byrow=TRUE) +
      (Y - P) / P / (1 - P)
    else {
      r <- X * matrix(cof, nrow=dx[1], ncol=dx[2], byrow=TRUE)
      R <- array(NA, dim=c(dx, k), dimnames=c(atx$dimnames, list(lev2)))
      for(j in 1:k) {
        y <- Y >= j
        p <- cumprob(L - cof.int[1] + cof.int[j])
        R[,,j] <- r + (y - p) / p / (1 - p)
      }
    }
    if(dopl) {
      xname <- atx$dimnames[[2]]; X <- unclass(X)
      for(i in 1:dx[2]) {
        if(pl == "loess") {
          if(k > 1)
            stop('pl="loess" not implemented for ordinal response')
          xi <- X[,i]
          ri <- r[,i]
          ddf <- data.frame(xi,ri)
          
          plot(xi, ri, xlim=if(missing(xlim)) range(xi) else xlim,
               ylim=if(missing(ylim)) range(ri) else ylim,
               xlab=xname[i], ylab=ylabpr)
          lines(lowess(xi,ri))
        }
        else if(k==1) {
          xi <- X[,i]; ri <- r[,i]
          plot(xi, ri, xlab=xname[i], ylab="Partial Residual",
               xlim=if(missing(xlim))range(xi) else xlim,
               ylim=if(missing(ylim))range(ri) else ylim)
          if(pl=="lowess") lines(lowess(xi, ri, iter=0, ...))
          else lines(supsmu(xi, ri, ...))
        }
        else {
          xi <- X[,i]
          ri <- R[,i,,drop=TRUE]
          smoothed <- vector('list',k)
          ymin <- 1e30; ymax <- -1e30
          
          for(j in 1:k) {
            w <- if(pl!='supsmu')
              lowess(xi, ri[,j], iter=0, ...)
            else
              supsmu(xi, ri[,j], ...)
            ymin <- min(ymin,w$y)
            ymax <- max(ymax,w$y)
            smoothed[[j]] <- w
          }
          plot(0, 0, xlab=xname[i], ylab=ylabpr,
               xlim=if(missing(xlim))range(xi) else xlim,
               ylim=if(missing(ylim))range(pretty(c(ymin,ymax)))
               else ylim, type='n')
          us <- par('usr')[1:2]
          for(j in 1:k) {
            w <- smoothed[[j]]
            lines(w, lty=j)
            if(is.character(label.curves)) {
              xcoord <- us[1]+(us[2]-us[1])*j/(k+1)
              text(xcoord, approx(w, xout=xcoord, rule=2, ties=mean)$y,
                   lev2[j])
            }
          }
          if(is.list(label.curves) || 
             (is.logical(label.curves) && label.curves))
            labcurve(smoothed, lev2, opts=label.curves)
        }
      }
      return(invisible(if(k==1)naresid(naa,r) else R))   
    }
    return(if(k==1) naresid(naa,r) else R)
  }
  
##if(type=='convexity') {
##  if(missing(p.convexity)) {
##	pq <- quantile(P, c(.01, .99))
##	if(pq[1]==pq[2]) pq <- range(P)
##	p.convexity <- seq(pq[1], pq[2], length=100)
## }
##  lcp <- length(p.convexity)
##  cp <- single(lcp)
##  for(i in 1:lcp) {
##	p <- p.convexity[i]
##	cp[i] <- mean(((p/P)^Y)*(((1-p)/(1-P))^(1-Y)))
##  }
##  if(pl) plot(p.convexity, cp, xlab='p', ylab='C(p)', type='l')
##  return(invisible(cp))
##}

  if(type %in% c("dfbeta", "dfbetas", "dffit", "dffits", "hat", "lp1")) {
    if(length(X <- unclass(object$x)) == 0)
      stop("you did not specify x=TRUE for the fit")
    v <- P * (1 - P)
    g <- lm(L + (Y - P) / v ~ X, weights=v)
    infl <- lm.influence(g)
    dfb <- coef(infl)    ## R already computed differences
    
    dimnames(dfb) <- list(rnam, c(cnam[kint], cnam[-(1:k)]))
    if(type=="dfbeta") return(naresid(naa, dfb))
    if(type=="dfbetas") {
      ## i <- c(kint, (k+1):length(cof))
      vv <- vcov(object, intercepts=1)
      return(naresid(naa, sweep(dfb, 2, diag(vv)^.5,"/")))
      ## was diag(object$var[i, i])
    }
    if(type=="hat") return(naresid(naa, infl$hat))
    if(type=="dffit") return(naresid(naa,
         (infl$hat * g$residuals)/(1 - infl$hat)))
    if(type=="dffits") return(naresid(naa,
         (infl$hat^.5)*g$residuals/(infl$sigma * (1 - infl$hat)) ))
    if(type=="lp1") return(naresid(naa,
         L - (infl$hat * g$residuals) / (1 - infl$hat)))
  }
}

residuals.orm <-
  function(object, 
           type=c("li.shepherd", "ordinary","score","score.binary","pearson",
             "deviance","pseudo.dep","partial",
             "dfbeta","dfbetas","dffit","dffits","hat","gof","lp1"),
           pl=FALSE, xlim, ylim, kint, label.curves=TRUE, 
           which, ...)
{
  type <- match.arg(type)
  args <- list(object=object, type=type, pl=pl,
               label.curves=label.curves, ...)
  if(!missing(kint))  args$kint  <- kint
  if(!missing(xlim))  args$xlim  <- xlim
  if(!missing(ylim))  args$ylim  <- ylim
  if(!missing(which)) args$which <- which
  do.call('residuals.lrm', args)
}

plot.lrm.partial <- function(..., labels, center=FALSE, ylim)
{
  dotlist <- list(...)
  nfit <- length(dotlist)
  if(missing(labels)) labels <- (as.character(sys.call())[-1])[1:nfit]

  vname <- dimnames(dotlist[[1]]$x)[[2]]
  nv <- length(vname)
  if(nv==0) stop('you did not specify x=TRUE on the fit')

  r <- vector('list', nv)
  for(i in 1:nfit) r[[i]] <- resid(dotlist[[i]], 'partial')

  for(i in 1:nv) {
    curves <- vector('list',nfit)
    ymin <- 1e10; ymax <- -1e10
    for(j in 1:nfit) {
      xx <- dotlist[[j]]$x[,vname[i]]
      yy <- r[[j]][,vname[i]]
      if(center)yy <- yy - mean(yy)
      curves[[j]] <- lowess(xx, yy, iter=0)
      ymin <- min(ymin, curves[[j]]$y)
      ymax <- max(ymax, curves[[j]]$y)
    }
    for(j in 1:nfit) {
      if(j==1) plot(curves[[1]], xlab=vname[i], ylab=NULL,
           ylim=if(missing(ylim)) c(ymin, ymax) else ylim, type='l')
      else lines(curves[[j]], lty=j)
    }
    if(nfit>1) labcurve(curves, labels)
  }
  invisible()
}
