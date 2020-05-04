psm <- function(formula=formula(data),
                data=parent.frame(),
                weights, subset, na.action=na.delete, dist='weibull', 
                init=NULL,  scale=0,
                control=survreg.control(),
                parms=NULL, model=FALSE, x=FALSE, y=TRUE, time.inc, ...) {
  
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  if(dist == 'extreme')
    warning('Unlike earlier versions of survreg, dist="extreme" does not fit\na Weibull distribution as it uses an identity link.  To fit the Weibull\ndistribution use the default for dist or specify dist="weibull".')
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
              names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action  ## FEH
  m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  special <- c("strata", "cluster")
  Terms <-
    if(missing(data)) terms(formula, special)
    else
      terms(formula, special, data=data)
  m$formula <- Terms
  ## Start FEH
  dul <- .Options$drop.unused.levels
  if(!length(dul) || dul) {
    on.exit(options(drop.unused.levels=dul))
    options(drop.unused.levels=FALSE)
  }

  m <- Design(eval.parent(m))
  atrx <- attributes(m)
  sformula <- atrx$sformula
  nact <- atrx$na.action
  Terms <- atrx$terms
  atr   <- atrx$Design
  ## End FEH
  
  weights <- model.extract(m, 'weights')
  Y <- model.extract(m, "response")
  Ysave <- Y
  
  ## Start FEH
  atY <- attributes(Y)
  
  ncy <- ncol(Y)
  maxtime <- max(Y[, - ncy])
  nnn <- c(nrow(Y), sum(Y[, ncy]))
  time.units <- units(Y)
  if(!length(time.units) || time.units == '') time.units <- "Day"
  if(missing(time.inc)) {
    time.inc <- switch(time.units, Day=30, Month=1, Year=1, maxtime / 10)
    if(time.inc >= maxtime | maxtime / time.inc > 25)
      time.inc <- max(pretty(c(0, maxtime))) / 10
  }
  ## End FEH

  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  strats <- attr(Terms, "specials")$strata
  cluster<- attr(Terms, "specials")$cluster
  dropx <- NULL
  if (length(cluster)) {
    if (missing(robust)) robust <- TRUE
    tempc <- untangle.specials(Terms, 'cluster', 1 : 10)
    ord <- attr(Terms, 'order')[tempc$terms]
    if (any(ord > 1)) stop ("Cluster can not be used in an interaction")
    cluster <- strata(m[, tempc$vars], shortlabel=TRUE)  #allow multiples
    dropx <- tempc$terms
  }
  if (length(strats)) {
    temp <- untangle.specials(Terms, 'strata', 1)
    dropx <- c(dropx, temp$terms)
    if (length(temp$vars) == 1) strata.keep <- m[[temp$vars]]
    else strata.keep <- strata(m[, temp$vars], shortlabel=TRUE)
    strata <- as.numeric(strata.keep)
    nstrata <- max(strata)
  }
  else {
    nstrata <- 1
    strata <- 0
  }
  
  if (length(dropx)) newTerms <- Terms[-dropx]
  else
    newTerms <- Terms
  X <- model.matrix(newTerms,m)
  
  ## Start FEH
  rnam <- dimnames(Y)[[1]]
  dimnames(X) <- list(rnam, c("(Intercept)", atr$colnames))
  ## End FEH except for 23nov02 and later changes
  
  n <- nrow(X)
  nvar <- ncol(X)
  
  offset <- atrx$offset
  if(!length(offset)) offset <- rep(0, n)
  
  if (is.character(dist)) {
    dlist <- survreg.distributions[[dist]]
    if (is.null(dlist)) stop(paste(dist, ": distribution not found"))
  }
  else
    if (is.list(dist)) dlist <- dist
    else
      stop("Invalid distribution object")
  if (!length(dlist$dist)) {
    ## WHAT IS THIS?
    if (is.character(dlist$name) && is.function(dlist$init) &&
        is.function(dlist$deviance)) {}
    else stop("Invalid distribution object")
  }
  else {
    if (!is.character(dlist$name) || !is.character(dlist$dist) ||
        !is.function(dlist$trans) || !is.function(dlist$dtrans))
      stop("Invalid distribution object")
  }	
  
  type <- attr(Y, "type")
  if (type ==  'counting') stop ("Invalid survival type")
  
  logcorrect <- 0   #correction to the loglik due to transformations
  if (length(dlist$trans)) {
    tranfun <- dlist$trans
    exactsurv <- Y[,ncol(Y)]  == 1
    if (any(exactsurv))
      logcorrect <-
        ifelse(length(weights),
               sum(weights[exactsurv]*logb(dlist$dtrans(Y[exactsurv, 1]))),
               sum(logb(dlist$dtrans(Y[exactsurv, 1]))))
    
    if (type == 'interval') {
      if (any(Y[,3] == 3))
        Y <- cbind(tranfun(Y[,1:2]), Y[,3])
      else Y <- cbind(tranfun(Y[,1]), Y[,3])
    }
    else {
      if (type == 'left')
        Y <- cbind(tranfun(Y[, 1]), 2 - Y[, 2])
      else
        Y <- cbind(tranfun(Y[, 1]), Y[, 2])
    }
    if (!all(is.finite(Y))) 
      stop("Invalid survival times for this distribution")
  }
  else {
    if (type == 'left') Y[, 2] <- 2- Y[, 2]
    else
      if (type == 'interval' && all(Y[, 3] < 3)) Y <- Y[, c(1, 3)]
  }
  
  ## if (!length(dlist$itrans)) itrans <- function(x) x
  ## else
  ##   itrans <- dlist$itrans
  
  if (length(dlist$scale)) {
    if (!missing(scale))
      warning(paste(dlist$name, 
                    "has a fixed scale, user specified value ignored"))
    scale <- dlist$scale
  }
  if (length(dlist$dist)) dlist <- survreg.distributions[[dlist$dist]]
  
  if (missing(control)) control <- survreg.control(...)
  
  if (scale < 0) stop("Invalid scale value")
  if (scale >0 && nstrata >1) 
    stop("Cannot have multiple strata with a fixed scale")
  
  ## Check for penalized terms
  pterms <- sapply(m, inherits, 'coxph.penalty')
  if (any(pterms)) {
    pattr <- lapply(m[pterms], attributes)
    ## 
    ## the 'order' attribute has the same components as 'term.labels'
    ##   pterms always has 1 more (response), sometimes 2 (offset)
    ## drop the extra parts from pterms
    temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
    if (length(dropx)) temp <- c(temp, dropx+1)
    pterms <- pterms[-temp]
    temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
    ord <- attr(Terms, 'order')[temp]
    if (any(ord > 1)) stop ('Penalty terms cannot be in an interaction')
    ##pcols <- (attr(X, 'assign')[-1])[pterms]
    assign <- attrassign(X,newTerms)
    pcols <- assign[-1][pterms]
    
    fit <- survpenal.fit(X, Y, weights, offset, init=init,
                         controlvals = control,
                         dist= dlist, scale=scale,
                         strata=strata, nstrat=nstrata,
                         pcols, pattr,assign, parms=parms)
  }
  else fit <-
    survreg.fit(X, Y, weights, offset, 
                init=init, controlvals=control,
                dist= dlist, scale=scale, nstrat=nstrata, 
                strata, parms=parms)
  
  if (is.character(fit))
    fit <- list(fail=fit)  #error message
  else {
    if (scale == 0) {
      nvar <- length(fit$coef) - nstrata
      fit$scale <- exp(fit$coef[-(1:nvar)])
      if (nstrata == 1) names(fit$scale) <- NULL
      else names(fit$scale) <- levels(strata.keep)
      fit$coefficients  <- fit$coefficients[1:nvar]
      fit$idf  <- 1 + nstrata
    }
    else {
      fit$scale <- scale
      fit$idf  <- 1
    }
    fit$loglik <- fit$loglik + logcorrect
  }
  
  if(length(nact)) fit$na.action <- nact  ## FEH
  fit$df.residual <- n - sum(fit$df)
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$means <- apply(X,2, mean)
  fit$call <- call
  fit$sformula <- sformula
  fit$dist <- dist
  fit$df.resid <- n-sum(fit$df) ##used for anova.survreg
  if (model) fit$model <- m
  if (x)     fit$x <- X[, -1, drop=FALSE]
  ##if (y)     fit$y <- Y   #FEH
  if (length(parms)) fit$parms <- parms
  ## Start FEH
  ##if (any(pterms)) class(fit)<- c('survreg.penal', 'survreg')
  ##else class(fit) <- 'survreg'
  fit$assign <- DesignAssign(atr, 1, Terms)
  fit$formula <- formula
  if(y) {
    class(Ysave) <- 'Surv'
    attr(Ysave, 'type') <- atY$type
    fit$y <- Ysave
  }
  scale.pred <-
    if(dist %in% c('weibull','exponential','lognormal','loglogistic'))
      c('log(T)','Survival Time Ratio')
    else 'T'
  
  logtest <- 2 * diff(fit$loglik)
  Nn <- if(length(weights)) sum(weights) else nnn[1]
  R2.max <- 1 - exp(2. * fit$loglik[1] / Nn)
  R2 <- (1 - exp(-logtest/Nn)) / R2.max
  df <- length(fit$coef) - 1
  P  <- if(df == 0) NA else 1. - pchisq(logtest, df)
  gindex <- GiniMd(fit$linear.predictors)
  Dxy <- if(type %in% c('right', 'left'))
           dxy.cens(fit$linear.predictors, Y)['Dxy']
  else {
    warning('Dxy not computed since right or left censoring not in effect')
    NA
  }
  stats <- c(nnn, logtest, df, P, R2, Dxy, gindex, exp(gindex))
  names(stats) <- c("Obs", "Events", "Model L.R.", "d.f.", "P",
                    "R2", "Dxy", "g", "gr")
  if(length(weights)) stats <- c(stats, 'Sum of Weights'=sum(weights))
  fit <- c(fit, list(stats=stats, weights=weights,
                     maxtime=maxtime, units=time.units,
                     time.inc=time.inc, scale.pred=scale.pred,
                     non.slopes=1, Design=atr, fail=FALSE))
  class(fit) <-
    if (any(pterms)) c('psm','rms','survreg.penal','survreg')
    else
      c('psm','rms','survreg')
  ## End FEH
  
  fit
  
}

Hazard   <- function(object, ...) UseMethod("Hazard")
Survival <- function(object, ...) UseMethod("Survival")

Hazard.psm <- function(object, ...) {
  dist <- object$dist
  g <- survreg.auxinfo[[dist]]$hazard
  formals(g) <- list(times=NA, lp=NULL, parms=logb(object$scale))
  g
}

Survival.psm <- function(object, ...) {
  dist <- object$dist
  g <- survreg.auxinfo[[dist]]$survival
  formals(g) <- list(times=NULL, lp=NULL, parms=logb(object$scale))
  g
}

Quantile.psm <- function(object, ...) {
  dist <- object$dist
  g <- survreg.auxinfo[[dist]]$Quantile
  formals(g) <- list(q=.5, lp=NULL, parms=logb(object$scale))
  g
}


Mean.psm <- function(object, ...) {
  dist <- object$dist
  g <- survreg.auxinfo[[dist]]$mean
  formals(g) <- list(lp=NULL, parms=logb(object$scale))
  g
}

predict.psm <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms","ccterms","adjto",
             "adjto.data.frame",  "model.frame"),
           se.fit=FALSE, conf.int=FALSE,
           conf.type=c('mean','individual','simultaneous'),
           kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...) {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               kint=1, na.action=na.action, expand.na=expand.na,
               center.terms=center.terms, ...)
  }


residuals.psm <-
  function(object,
           type=c("censored.normalized",
             "response", "deviance","dfbeta","dfbetas", 
             "working","ldcase","ldresp","ldshape", "matrix", "score"), ...) {
  type <- match.arg(type)
  if(type != 'censored.normalized') {
    r <- getS3method('residuals', 'survreg')
    s <- if(type == 'score') {
      X <- cbind('(Intercept)'=1, object$x)
      if(! length(X)) stop('did not use x=T with fit')
      wts <- object$weights
      if(! length(wts)) wts <- 1
      res <- r(object, type='matrix')
      s <- as.vector(res[, 'dg']) * wts * X
      if(NROW(object$var) > length(coef(object)))
        s <- cbind(s, 'Log(scale)'=unname(res[,'ds']))
    } else r(object, type=type)
    return(s)
  }
  
  y <- object$y
  aty <- attributes(y)
  if(length(y) == 0) stop('did not use y=T with fit')
  ncy <- ncol(y)
  scale <- object$scale
  dist  <- object$dist
  trans <- survreg.distributions[[dist]]$trans
  r <- (trans(y[, -ncy, drop=FALSE]) - object$linear.predictors) / scale
  label(r) <- 'Normalized Residual'
  ev <- y[, ncy]
  lab <- aty$inputAttributes$event$label
  if(length(lab)) label(ev) <- lab
  ## Moved the following line here from bottom
  r <- Surv(r, ev)
  if(length(object$na.action)) r <- naresid(object$na.action, r)
  attr(r,'dist') <- dist
  attr(r,'type') <- aty$type
  class(r) <- c('residuals.psm.censored.normalized','Surv')
  g <- survreg.auxinfo[[dist]]$survival
  formals(g) <- list(times=NULL, lp=0, parms=0)
  attr(r,'theoretical') <- g
  r
}

lines.residuals.psm.censored.normalized <- 
  function(x, n=100, lty=1, xlim=range(r[,-ncol(r)],na.rm=TRUE),
           lwd=3, ...) {
  r <- x
  x <- seq(xlim[1], xlim[2], length=n)
  tx <- x
  dist <- attr(r, 'dist')
  if(dist %in% c('weibull','loglogistic','lognormal')) tx <- exp(x)
  ## $survival functions log x
  lines(x, attr(r,'theoretical')(tx), lwd=lwd, lty=lty, ...)
  invisible()
}

survplot.residuals.psm.censored.normalized <- 
  function(fit, x, g=4, col, main, ...) {
    r <- fit
    
    if(missing(x)) {
      survplot(npsurv(r ~ 1), conf='none', xlab='Residual', 
               col=if(missing(col))par('col') else col, ...)
      if(!missing(main)) title(main)
    }
  else {
    if(is.character(x)) x <- as.factor(x)
    if(!is.factor(x) && length(unique(x))>5) x <- cut2(x, g=g)
    s <- is.na(r[,1]) | is.na(x)
    if(any(s)) {r <- r[!s,]; x <- x[!s,drop=TRUE]}
    survplot(npsurv(r ~ x, data=data.frame(x,r)),  xlab='Residual',
             conf='none',
             col=if(missing(col))1:length(levels(x)) else par('col'), ...)
    if(missing(main)) {
      main <-
        if(length(lab <- attr(x,'label'))) lab 
        else ''
    }
    if(main != '') title(main)
  }
    lines(r, lty=1, lwd=3)
    invisible()
  }
