## ?? Why does confint call use nrp???
## Value adjusted to is irrelevant when the factor does not interact with
# other factors.  Form of factors is as follows: factor1=value1,factor2=val2:
# Values:
#	NA	: test factor, use all default settings
#	w	: adjust this factor to w when estimating effects of others
#	c(lo,hi): use range for effect (lo,hi), adjust to default value
#	c(lo,w,hi): use range (lo,hi), adjust to w.  Any of 3 can be NA.
# For categories and strata values can be character
# values that are original values before translating to factors -
# only enough letters are needed to uniquely identify the category 
# This applies to category and strata vars.  Default adjusted to is
# from second element of limits vector.
# For category factors, all comparisons to reference category are made.
# Reference category is assumed to be adjusted to value.
# est.all is T to estimate effects for all factors, not just those listed
# in ...

summary.rms <- function(object, ..., est.all=TRUE, antilog, conf.int=.95,
                        abbrev=FALSE, vnames=c("names","labels"),
                        conf.type=c('individual','simultaneous'),
                        usebootcoef=TRUE,
                        boot.type=c('percentile','bca','basic'),
                        verbose=FALSE)
{	
  obj.name <- as.character(sys.call())[2]
  at <- object$Design
  labels <- at$label

  vnames <- match.arg(vnames)
  conf.type <- match.arg(conf.type)
  boot.type <- match.arg(boot.type)
  blabel <- switch(boot.type,
                   percentile = 'bootstrap nonparametric percentile',
                   bca = 'bootstrap BCa',
                   basic = 'basic bootstrap')
  ## if(conf.type == 'simultaneous') require(multcomp)
  
  assume <- at$assume.code
  if(is.null(assume)) stop("fit does not have design information")
  if(any(assume==10))
    warning("summary.rms does not currently work with matrix factors in model")
  name  <- at$name
  parms <- at$parms

  scale <- object$scale.pred
  if(missing(antilog)) antilog <- length(scale)==2
  if(antilog & length(scale) < 2) scale <- c("","Antilog")
  ## Hazard ratio is antilog of negative of difference if orm and
  ## ratio pertains
  logRatioAdj <- if(inherits(object, 'orm') && object$family == 'loglog') -1 else 1

  factors <- rmsArgs(substitute(list(...)))
  nf <- length(factors)

  if(est.all) which <- (1:length(assume))[assume!=9]
  if(nf > 0) {
    jw <- charmatch(names(factors), name, 0)
    if(any(jw == 0)) stop(paste("factor name(s) not in the design:",
             paste(names(factors)[jw == 0], collapse=" ")))
    if(!est.all) which <- jw
    if(any(assume[which] == 9))
      stop("cannot estimate effects for interaction terms alone")
  }
  
  Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  values <- Limval$values
  ## The next statement (9Jun98) makes limits[1:3,] keep all levels of
  ## factors.  Problem is that [.data.frame does not pass drop to []
  ## when first subscripts are specified
  oldopt <- options(drop.factor.levels=FALSE)
  on.exit(options(oldopt))

  lims <- Limval$limits[1 : 3 , , drop=FALSE]
  
  ##Find underlying categorical variables
  ucat <- rep(FALSE, length(assume))
  for(i in (1:length(assume))[assume != 5 & assume < 8])
    ucat[i] <- name[i] %in% names(values) &&
  length(V <- values[[name[i]]]) && is.character(V)
  
  stats <- lab <- NULL
  lc <- length(object$coef)
  ##Number of non-slopes:
  nrp <- num.intercepts(object, 'coef')
  nrp1 <- nrp + 1
  ## Exclude non slopes
  beta <- object$coef[nrp1 : lc]
  var <- vcov(object, regcoef.only=TRUE, intercepts='none')

  zcrit <- if(length(idf <- object$df.residual)) qt((1 + conf.int)/2, idf)
   else qnorm((1 + conf.int) / 2)
  cll <- paste(signif(conf.int, 3))
  bcoef <- if(usebootcoef) object$boot.Coef
  if(length(bcoef)) bcoef <- bcoef[, nrp1 : lc, drop=FALSE]

  jf <- 0
  if(nf > 0) for(i in jw) {
    jf <- jf + 1
    z <- value.chk(at, i, factors[[jf]], 0, Limval)
    lz <- length(z)
    if(lz == 1 && !is.na(z)) lims[2, i] <-  z
    if(lz == 2) {
      if(!is.na(z[1])) lims[1, i] <- z[1]
      if(!is.na(z[2])) lims[3, i] <- z[2]
    }
    else
      if(lz == 3) lims[!is.na(z), i] <- z[!is.na(z)]
    if(lz < 1 | lz > 3) stop("must specify 1,2, or 3 values for a factor")
  }
  adj <- lims[2,, drop=FALSE]
  isna <- sapply(adj, is.na)

  if(any(isna))
    stop(paste("adjustment values not defined here or with datadist for",
               paste(name[assume != 9][isna], collapse=" ")))
  k <- which[assume[which] %nin% c(8, 5, 10) & !ucat[which]]
  m <- length(k)
  if(m) {
    isna <- is.na(lims[1, name[k], drop=FALSE] + lims[3, name[k], drop=FALSE])
    ##note char. excluded from k
    if(any(isna)) stop(paste("ranges not defined here or with datadist for",
                             paste(name[k[isna]], collapse=" ")))
  }
  
  xadj <- unclass(rms.levels(adj, at))
  if(m) {
    adj <- xadj
    M <- 2 * m
    odd  <- seq(1, M, by=2)
    even <- seq(2, M, by=2)
    ##Extend data frame
    for(i in 1:length(adj)) adj[[i]] <- rep(adj[[i]], M)
    
    i <- 0
    for(l in k) {
      i <- i + 1
      adj[[name[l]]][(2 * i - 1) : (2 * i)] <- lims[c(1, 3), name[l]]
    }
    xx <- predictrms(object, newdata=adj, type="x")
    xd <- matrix(xx[even,] - xx[odd,], nrow=m)
    xb <- xd %*% beta
    se <- drop((((xd %*% var) * xd) %*% rep(1, ncol(xd)))^.5)
    if(conf.type == 'simultaneous' && length(xb) > 1) {
      if(verbose) {
        cat('Confidence intervals are simultaneous for these estimates:\n')
        print(as.vector(xb))
      }
      u <- confint(multcomp::glht(object,
                        cbind(matrix(0, nrow=nrow(xd), ncol=nrp), xd),
                        df=if(length(idf)) idf else 0),
                   level=conf.int)$confint
      low <- u[, 'lwr']
      up  <- u[, 'upr']
    } else if(length(bcoef)) {
      best <- t(xd %*% t(bcoef))
      lim <- bootBCa(xb, best, type=boot.type, n=nobs(object), seed=object$seed,
                     conf.int=conf.int)
      if(is.matrix(lim)) {
        low <- lim[1,]
        up  <- lim[2,]
      } else {
        low <- lim[1]
        up <- lim[2]
      }
    } else {
      low <- xb - zcrit*se
      up  <- xb + zcrit*se
    }
    lm <- as.matrix(lims[, name[k], drop=FALSE])
    stats <- cbind(lm[1,], lm[3,], lm[3,] - lm[1,], xb, se, low, up, 1)
    lab <- if(vnames=='names') name[k] else labels[k]
    if(antilog) {
      stats <- rbind(stats,
                     cbind(stats[, 1 : 3,drop=FALSE],
                           exp(logRatioAdj * xb), NA, exp(low), exp(up), 2))
      lab <- c(lab, rep(paste("", scale[2]), m))
      w <- integer(M)
      w[odd] <- 1 : m
      w[even]<- m + (1 : m)
      stats <- stats[w,]
      lab <- lab[w]
    }
  }
  
  for(j in 1:length(xadj)) xadj[[j]] <- rep(xadj[[j]], 2)
  
  for(i in which[assume[which] == 5 | ucat[which]]) {
    ## All comparisons with reference category
    
    parmi <- if(ucat[i]) values[[name[i]]] else parms[[name[i]]]
    parmi.a <- if(abbrev) abbreviate(parmi) else parmi
    iref <- as.character(xadj[[name[i]]][1])
    ki <- match(iref, parmi)
    for(j in parmi) {
      if(j != iref) {
        kj <- match(j, parmi)
        adj <- xadj
        adj[[name[i]]] <- c(iref, j)
        adj <- as.data.frame(adj)
        xx <- predictrms(object, newdata=adj, type="x")
        xd <- matrix(xx[2,] - xx[1,], nrow=1)
        xb <- xd %*% beta
        se <- sqrt((xd %*% var) %*% t(xd))
        if(conf.type == 'simultaneous' && length(xb) > 1) {
          if(verbose) {
            cat('Confidence intervals are simultaneous for these estimates:\n')
            print(as.vector(xb))
          }
          u <- confint(multcomp::glht(object,
                            cbind(matrix(0, nrow=nrow(xd), ncol=nrp), xd),
                            df=if(length(idf)) idf else 0),
                       level=conf.int)$confint
          low <- u[,'lwr']
          up  <- u[,'upr']
        } else if(length(bcoef)) {
          best <- t(xd %*% t(bcoef))
          lim <- bootBCa(xb, best, type=boot.type,
                         n=nobs(object), seed=object$seed, conf.int=conf.int)
          if(is.matrix(lim)) {
            low <- lim[1,]
            up  <- lim[2,]
          } else {
            low <- lim[1]
            up  <- lim[2]
          }
        } else {
          low <- xb - zcrit*se
          up  <- xb + zcrit*se
        }
        stats <- rbind(stats,cbind(ki, kj, NA,
                                   xb, se, low, up, 1))
        lab <-c(lab,
                paste(if(vnames=='names') name[i] else labels[i],
                      " - ", parmi.a[kj], ":",
                      parmi.a[ki], sep=""))
        if(antilog) {
          stats <- rbind(stats,cbind(ki, kj, NA,
                                     exp(logRatioAdj * xb),
                                     NA, exp(low), exp(up), 2))
          lab <- c(lab, paste("", scale[2]))}
      }
    }
  }
  
  dimnames(stats) <-
    list(lab, c("Low", "High",
                "Diff.", "Effect", "S.E.",
                paste("Lower", cll), paste("Upper", cll), "Type"))
  
  attr(stats, "heading") <-
    paste("             Effects              Response : ",
          as.character(formula(object))[2], sep='')
  attr(stats,"class") <- c("summary.rms", "matrix")
  attr(stats,"scale") <- scale
  attr(stats,"obj.name") <- obj.name
  interact <- at$interactions
  adjust <- ""
  if(length(interact)) {
    interact <- sort(unique(interact[interact > 0]))
    nam <- name[which[match(which, interact, 0) > 0]]
    if(length(nam)) for(nm in nam) 
      adjust <- paste(adjust, nm, "=",
                      if(is.factor(xadj[[nm]]))
                      as.character(xadj[[nm]])[1] else
                      format(xadj[[nm]][1]), " ", sep="")
  }
  attr(stats, "adjust") <- adjust
  attr(stats, "conf.type") <-
    if(length(bcoef)) blabel else 'z'
  
  stats
}

print.summary.rms <- function(x, ...)
{
  cstats <- dimnames(x)[[1]]
  for(i in 1 : 7) cstats <- cbind(cstats, format(signif(x[, i], 5)))
#  for(i in 4 : 7) cstats <- cbind(cstats, format(round (x[, i], 2)))
  dimnames(cstats) <- list(rep("", nrow(cstats)), 
                           c("Factor", dimnames(x)[[2]][1 : 7]))
  cat(attr(x,"heading"), "\n\n")
  print(cstats, quote=FALSE)
  if((A <- attr(x, "adjust")) != "") cat("\nAdjusted to:", A, "\n")
  blab <- switch(attr(x, 'conf.type'),
                 'bootstrap nonparametric percentile' = 
                  'Bootstrap nonparametric percentile confidence intervals',
                 'bootstrap BCa' = 'Bootstrap BCa confidence intervals',
                 'basic bootstrap' = 'Basic bootstrap confidence intervals',
                 '')
  if(blab != '') cat('\n', blab, '\n', sep='')
  cat('\n')
  invisible()
}


latex.summary.rms <-
  function(object, 
           title=paste('summary', attr(object, 'obj.name'), sep='.'),
           table.env=TRUE, ...)
{ 

  title <- title   # because of lazy evaluation
  caption <- latexTranslate(attr(object, "heading"))
  scale <- attr(object, "scale")
  object <- object[, -8, drop=FALSE]
  rowl <- latexTranslate(dimnames(object)[[1]])
  rowl <- ifelse(substring(rowl, 1, 1) == " ",
                 paste("~~{\\it ",
                       substring(rowl,2), "}", sep=""),
                 rowl) # preserve leading blank
  rowl <- sedit(rowl, "-", "---")
  cstats <- matrix("", nrow=nrow(object), ncol=ncol(object), 
                   dimnames=dimnames(object))
  for(i in 1 : 7) cstats[,i] <- format(signif(object[, i], 5))
  ## for(i in 4 : 7) cstats[,i] <- format(round(object[, i],  2))
  cstats[is.na(object)] <- ""
  caption <- sedit(caption, "    Response","~~~~~~Response")
  cstats <- as.data.frame(cstats)
  attr(cstats,"row.names") <- rowl
  names(cstats)[3] <- "$\\Delta$"
  latex(cstats, caption=if(table.env) caption else NULL,
        title=title, rowlabel="",
        col.just=rep("r", 7), table.env=table.env, ...)
}


# was q=c(.7, .8, .9, .95, .99)
plot.summary.rms <-
  function(x, at, log=FALSE, 
           q=c(0.9, 0.95, 0.99), xlim, nbar, cex=1, nint=10,
           cex.main=1, clip=c(-1e30, 1e30), main,
           col=rgb(red=.1, green=.1, blue=.8, alpha=c(.1,.4,.7)),
           col.points=rgb(red=.1, green=.1, blue=.8, alpha=1),
           pch=17, lwd=if(length(q) == 1) 3 else 2 : (length(q) + 1), ...)
{
  confbar <-
    function(y, est, se, q, col, col.points,
             pch=17, lwd=rep(3, length(q)), clip=c(-1e30, 1e30),
             fun  = function(x) x, 
             qfun = function(x) ifelse(x==.5, qnorm(x),
               ifelse(x < .5, qnorm(x / 2), qnorm((1 + x) / 2)))) {

      n <- length(q)
      q <- c(1 - rev(q), .5, q)
      a <- fun(est)
      points(a, y, col=col.points, pch=pch)
      a <- fun(est + se * qfun(q))
      a[a < clip[1]] <- NA; a[a > clip[2]] <- NA
      m <- length(q)
      segments(c(a[1], a[m]), y, c(a[2], a[m - 1]), y, col=col[1], lwd=lwd[1])
      if(n > 1) segments(c(a[2], a[m - 1]), y, c(a[3], a[m - 2]),
                         col=col[2], lwd=lwd[2])
      if(n > 2) segments(c(a[3], a[m - 2]), y, c(a[4], a[m - 3]),
                         col=col[3], lwd=lwd[3])
      names(a) <- format(q)
      invisible(a)
    }

  scale  <- attr(x, "scale")
  adjust <- attr(x, "adjust")

  Type   <- x[, "Type"]
  x  <- x[Type==1,, drop=FALSE]
  lab    <- dimnames(x)[[1]]
  effect <- x[, "Effect"]
  se     <- x[, "S.E."]
  if(!log && any(Type==2)) {
    fun <- exp
    tlab <- scale[2]
  }
  else {
    fun <- function(x) x
    if(log) {
      if(length(scale) == 2) tlab <- scale[2]
      else tlab <- paste("exp(", scale[1], ")", sep="")
    }
    else tlab <- scale[1]
  }
  if(!length(scale)) tlab <- ''  ## mainly for Glm fits
  if(!missing(main)) tlab <- main
  augment <- if(log | any(Type==2)) c(.1, .5, .75, 1) else 0
  n     <- length(effect)
  out   <- qnorm((max(q) + 1) / 2)
  if(missing(xlim) && !missing(at))
    xlim <- range(if(log) logb(at) else at) else
  if(missing(xlim)) {
    xlim <- fun(range(c(effect - out * se, effect + out * se)))
    xlim[1] <- max(xlim[1], clip[1])
    xlim[2] <- min(xlim[2], clip[2])
  }
  else
    augment <- c(augment, if(log) exp(xlim) else xlim)
  
  fmt <- function(k) {
    m <- length(k)
    f <- character(m)
    for(i in 1 : m) f[i] <- format(k[i])
    f
  }
  lb <- ifelse(is.na(x[, 'Diff.']), lab,
               paste(lab, ' - ',
                     fmt(x[, 'High']), ':', fmt(x[, 'Low']), sep=''))
  plot.new(); par(new=TRUE)
  mxlb <- .1 + max(strwidth(lb, units='inches', cex=cex))
  tmai <- par('mai')
  on.exit(par(mai=tmai))
  par(mai=c(tmai[1], mxlb, 1.5*tmai[3], tmai[4]))
  
  outer.widths <- fun(effect + out * se)-fun(effect - out * se)
  if(missing(nbar)) nbar <- n
  npage <- ceiling(n/nbar)
  is <- 1
  for(p in 1 : npage) {
    ie <- min(is + nbar - 1, n)
    plot(1:nbar, rep(0,nbar), xlim=xlim, ylim=c(1,nbar),
         type="n", axes=FALSE, 
         xlab="", ylab="")
    if(cex.main > 0) title(tlab, cex=cex.main)
    lines(fun(c(0, 0)), c(nbar - (ie - is), nbar), lty=2)
    if(log) {
      pxlim <- pretty(exp(xlim), n=nint)
      pxlim <- sort(unique(c(pxlim, augment)))
      ## For wome weird reason, sometimes duplicates (at xlim[2])
      ## still remain
      pxlim <- pxlim[pxlim >= exp(xlim[1])]
      if(!missing(at)) pxlim <- at
      axis(3, logb(pxlim), labels=format(pxlim))
    }
    else {
      pxlim <- pretty(xlim, n=nint)
      pxlim <- sort(unique(c(pxlim, augment)))
      pxlim <- pxlim[pxlim >= xlim[1]]
      if(!missing(at)) pxlim <- at
      axis(3, pxlim)
    }
    imax <- (is : ie)[outer.widths[is : ie] == max(outer.widths[is : ie])][1]
    for(i in is : ie) {
      confbar(nbar - (i - is + 1) + 1, effect[i], se[i], q=q,
              col=col, col.points=col.points,
              fun=fun, clip=clip, lwd=lwd, pch=pch)
      mtext(lb[i], 2, 0, at=nbar - (i - is + 1) + 1, cex=cex,
            adj=1, las=1)
    }
    if(adjust != "") {
      adjto <- paste("Adjusted to:", adjust, sep="")
      xx <- par('usr')[2]
      if(nbar > ie) text(xx, nbar - (ie - is + 1), adjto, adj=1, cex=cex)
      else title(sub=adjto, adj=1, cex=cex)
    }
    is <- ie + 1
  }
  invisible()
}
