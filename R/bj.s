bj <- function(formula=formula(data), data,
               subset, na.action=na.delete, 
               link="log",	control=NULL,
               method='fit', x=FALSE, y=FALSE, time.inc)
{
  call <- match.call()
  m <- match.call(expand.dots=FALSE)
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
              names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  dul <- .Options$drop.unused.levels
  if(!length(dul) || dul)
    {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }
  
  X <- Design(eval.parent(m))
  if(method=='model.frame') return(X)
  atrx  <- attributes(X)
  nact  <- atrx$na.action
  Terms <- atrx$terms
  atr   <- atrx$Design
  sformula <- atrx$sformula
  
  lnames <- c("logit","probit","cloglog","identity","log","sqrt",
              "1/mu^2","inverse")

  link <- pmatch(link, lnames, 0)
  if(link==0) stop("invalid link function")
  link <- lnames[link]
  Y <- model.extract(X, "response")
  atY <- attributes(Y)
  ncy <- ncol(Y)
  maxtime <- max(Y[,-ncy])
  nnn <- c(nrow(Y),sum(Y[,ncy]))
  if (!inherits(Y, "Surv")) stop("Response must be a survival object")
  
  type <- attr(Y, "type")
  
  linkfun <- make.link(link)$linkfun 
  
  if (type != 'right') stop ("Surv type must by 'right' censored")
  Y <- cbind(linkfun(Y[,1]), Y[,2])
  
  X <- model.matrix(sformula, X)
  assgn <- DesignAssign(atr, 1, Terms)
  
  if(method=='model.matrix') return(X)
  
  time.units <- units(Y)
  if(is.null(time.units) || time.units=='') time.units <- "Day"
  if(missing(time.inc))  {
    time.inc <- switch(time.units,Day=30,Month=1,Year=1,maxtime/10)
    if(time.inc >= maxtime | maxtime/time.inc > 25) 
      time.inc <- max(pretty(c(0, maxtime)))/10
  }
  rnam <- dimnames(Y)[[1]]
  dimnames(X) <- list(rnam, c("(Intercept)",atr$colnames))
  
  n <- nrow(X)
  nvar <- ncol(X)
  
  fit <- bj.fit(X, Y, control=control)
  if(link == 'log') fit$stats <- c(fit$stats, gr=unname(exp(fit$stats['g'])))
  
  if(fit$fail) {
    cat("Failure in bj.fit\n")
    return(fit)
  }
  
  if (length(nact)) fit$na.action <- nact
  
  fit <- c(fit, list(maxtime=maxtime, units=time.units,
                     time.inc=time.inc, non.slopes=1, assign=assgn))
  class(fit) <-  c("bj", "rms")
  fit$sformula <- sformula
  fit$terms   <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  fit$call    <- call
  fit$Design  <- atr
  if (x) fit$x <- X[, -1, drop=FALSE]
  if (y)
    {
      class(Y) <- 'Surv'
      attr(Y,'type') <- atY$type
      fit$y <- Y
    }
  scale.pred <- if(link=="log") c("log(T)","Survival Time Ratio") else "T"
  fit$scale.pred <- scale.pred
  fit$link       <- link
  fit
}

bj.fit <- function(x, y, control = NULL) {
  if(ncol(y) != 2)
	stop("y is not a right-censored Surv object")
  status <- y[, 2]
  yy <- y[, 1]
  iter.max <- control$iter.max
  eps <- control$eps
  trace <- control$trace
  tol <- control$tol
  max.cycle <- control$max.cycle
  if(length(iter.max) == 0) iter.max <- 20
  if(length(eps) == 0)   eps <- 0.001
  if(length(trace) == 0) trace <- FALSE
  if(length(tol) == 0)   tol <- 1e-007
  if(length(max.cycle) == 0) max.cycle <- 30
  x <- as.matrix(x)
  if(all(x[, 1] == 1)) x <- x[, -1, drop = FALSE]
  d <- dim(x)
  nvar <- d[2]
  if(length(nvar) == 0)	nvar <- 0
  N <- length(yy)
  if(nvar > 0)
    {
      xbar <- apply(x, 2, mean)
      xm <- x - rep(xbar, rep(N, nvar))
    }
  else xm <- 0
  timeorig <- yy
  order.orig <- 1:N
  dummystrat <- factor(rep(1, N))
  betahat <- rep(0, max(nvar, 1))
  betamatrix <- NULL
  sse <- 0
  n <- 0	##
  ## new stuff
  nonconv <- FALSE	##
  repeat
    {
      oldbeta <- betahat
      oldsse <- sse
      if(nvar == 0) ypred <- 0
	  else
        {
          betahat <- solvet(t(xm) %*% xm, t(xm) %*% yy, tol = tol)
          ypred <- x %*% betahat
        }
      alphahat <- mean(yy - ypred)
      sse <- sum((yy - ypred)^2)
      razlika <- oldsse/sse
      if(trace)
        cat("iteration = ", n, "   sse ratio = ", format(razlika), "\n")
      n <- n + 1
      if(trace)
        cat("  alpha = ", format(alphahat),
            "   beta = ", format(betahat), "\n\n")
      ehat <- timeorig - ypred
      if(!nonconv)
        {
          if(abs(razlika - 1) <= eps)
            break
          else if(n > iter.max)
            {
              cyclesse <- NULL
              cycleperiod <- 0
              nonconv <- TRUE
              firstsse <- sse
            }
        }
	  else
        {
          betamatrix <- cbind(betamatrix, c(alphahat, betahat))
          cyclesse <- c(cyclesse, sse)
          cycleperiod <- cycleperiod + 1
          if(any(abs(firstsse - cyclesse) < 1e-007))
            {
              cat("\nCycle period = ", cycleperiod, "\n")
              meanbeta <- apply(betamatrix, 1, mean)
              alphahat <- meanbeta[1]
              betahat <- meanbeta[2:length(meanbeta)]
              ypred <- x %*% betahat
              ehat <- timeorig - ypred
              break
            }
		  else if(cycleperiod >= max.cycle)
			break
        }
      state <- status
      state[ehat == max(ehat)] <- 1
      S <- structure(cbind(ehat, state), class = "Surv", type = "right")
      KM.ehat <-
        survfitKM(dummystrat, S, conf.type = "none", se.fit = FALSE)
      n.risk <- KM.ehat$n.risk
      surv <- KM.ehat$surv
      repeats <- c(diff( - n.risk), n.risk[length(n.risk)])
      surv <- rep(surv, repeats)
      w <-  - diff(c(1, surv))
      m <- order(ehat,  - status)
      bla <- cumsum((w * ehat[m]))
      bla <- (bla[length(bla)] - bla)/(surv + state[m])
      ## Put bla back into original order
      bl <- bla
      bl[(1:N)[m]] <- bla
      yhat <- if(nvar == 0) bl else x %*% betahat + bl
      yy[state == 0] <- yhat[state == 0]
    }
  n <- n - 1
  if(nonconv)
    {
      if(cycleperiod < max.cycle)
        cat("\nNo convergence in", n, "steps, but cycle found - average beta returned\n")
	  else
        {
          cat("\nNo convergence in", n, "steps\n")
          return(list(fail = TRUE))
        }
    }
  f <- list(fail = FALSE, iter = n)
  cof <- if(nvar == 0) alphahat else c(alphahat, betahat)
  dx <- dimnames(x)[[2]]
  if(length(dx) == 0 && nvar > 0) dx <- paste("x", 1:nvar, sep = "")
  names(cof) <- c("Intercept", dx)
  f$coefficients <- cof
  ehat.u <- ehat[status == 1]
  edf <- sum(status) - nvar - 1
  sigma <- sqrt(sum((ehat.u - mean(ehat.u))^2)/edf)
  if(nvar > 0)
    {
      x <- cbind(Intercept = 1, x)[status == 1,  , drop = FALSE]
      f$var <- solvet(t(x) %*% x, tol = tol) * sigma * sigma
    }
  else f$var <- (sigma * sigma)/N
  f$linear.predictors <- alphahat + as.vector(ypred)
  g <- GiniMd(f$linear.predictors)
  stats <- c(N, sum(status), nvar, edf, sigma, g)
  names(stats) <- c("Obs", "Events", "d.f.", "error d.f.", "sigma", "g")
  f$stats <- stats
  if(any(status == 0))
	yy <- structure(yy, class = "impute", imputed = (1:N)[status == 0])
  f$y.imputed <- yy
  f
}

bjplot <- function(fit, which=1:dim(X)[[2]])
{
  if(!all(c('x','y') %in% names(fit)))
	stop('must specify x=TRUE,y=TRUE to bj to use bjplot')
  X <- fit$x
  Y <- fit$y
  xnam <- dimnames(X)[[2]]
  yy <- fit$y.imputed
  imp <- is.imputed(yy)
  trans <- if(fit$link=='identity') '' else fit$link

  ## Do Hillis plot first
  N <- length(fit$y[, 1])
  dummystrat <- factor(rep(1, N))
  S <- resid(fit)
  S[S[, 1] == max(S[, 1]), 2] <- 1
  m <- order(fit$y[, 1],  - fit$y[, 2])
  resd <- S[m, 1]
  cens <- S[m, 2]
  KM.ehat <- survfitKM(dummystrat, S, 
                       conf.type = "none", se.fit = FALSE)
  repeats <- c(diff( - KM.ehat$n.risk), KM.ehat$n.risk[length(KM.ehat$n.risk)])
  if(length(KM.ehat$time) != N)
    {
      time <- rep(KM.ehat$time, repeats)
      surv <- rep(KM.ehat$surv, repeats)
    }
  else
    {
	  time <- KM.ehat$time
	  surv <- KM.ehat$surv
	}
  u <- runif(N-1, 0, surv[1:(N - 1)])
  w <- approx(surv, time, xout=u, method='constant', f=0)
  t.i <- c(w$y, max(time))
  surv.i <- c(w$x, min(surv))
  residnew <- resd
  residnew[cens == 0] <- t.i[cens == 0]
  retlist <- list(predictor = fit$linear.predictor[m], 
				  x = fit$x[m,  ], res.cens = resd, hillis = residnew, 
				  cens = cens)
  predictor <- fit$linear.predictor[m]
  plot(predictor, resd, type = "n", 
	   xlab = "Linear Predictor", ylab = "Residuals")
  points(predictor[cens == 0], resd[cens == 0], pch = 1)
  points(predictor[cens == 1], resd[cens == 1], pch = 16)
  plot(predictor, residnew, type = 	"n", xlab = "Linear Predictor", 
	   ylab = "Residuals")
  points(predictor[cens == 0], residnew[cens == 0], pch = 1)
  points(predictor[cens == 1], residnew[cens == 1], pch = 16)
  
  for(i in which)
    {
      xi <- X[,i]
      ry <- range(yy,Y)

      plot(xi, Y[,1], xlab=xnam[i], ylab=paste('Observed',trans,'Time'),
           type='n', ylim=ry)
      points(xi[!imp], Y[!imp,1], pch=16)
      if(any(imp))
        {
          points(xi[imp],  Y[imp,1],  pch=1)

          plot(xi, yy, xlab = xnam[i], ylab=paste('Imputed',trans,'Time'), 
               type = "n", ylim=ry)
          points(xi[imp],   yy[imp],  pch = 1)
          segments(xi[imp], Y[imp,1], xi[imp], yy[imp])
          points(xi[!imp],  yy[!imp], pch = 16)

          plot(xi, yy, xlab=xnam[i],
               ylab=paste('Observed or Imputed',trans,'Time'),
               type='n', ylim=ry)
          points(xi[!imp], yy[!imp], pch=16)
          points(xi[imp],  yy[imp],  pch=1)
        }
    }
  invisible(retlist)
}

print.bj <- function(x, digits=4, long=FALSE, coefs=TRUE, latex=FALSE,
                     title="Buckley-James Censored Data Regression", ...)
{
  k <- 0
  z <- list()
  

  if(length(zz <- x$na.action))
    {
      k <- k + 1
      z[[k]] <- list(type=paste('naprint', class(zz)[1], sep='.'), list(zz))
    }
  
  stats <- x$stats
  ci    <- x$clusterInfo
  misc   <- reVector(Obs          = stats['Obs'],
                     Events       = stats['Events'],
                     'Cluster on' = ci$name,
                     'Clusters'   = ci$n)
  dfstat <- reVector('Regression d.f.' = stats['d.f.'],
                     sigma=stats['sigma'],
                     'd.f.'=stats['error d.f.'])
  disc <- reVector(g = stats['g'], gr = stats['gr'])
  k <- k + 1
  z[[k]] <- list(type='stats',
                 list(headings=list('', '', c('Discrimination','Indexes')),
                      data=list(misc, c(dfstat,c(NA,digits,NA)), c(disc, 3))))
  
  cof <- x$coefficients
  se <- sqrt(diag(x$var))
  k <- k + 1
  z[[k]] <- list(type='coefmatrix',
                 list(coef = cof, se = se))

  p <- length(cof)
  if(long &&  p > 1)
    {
	  ss <- diag(1/se)
	  correl <- ss %*% x$var %*% ss
	  dimnames(correl) <- list(names(cof), names(cof))
	  ll <- lower.tri(correl)
	  correl[ll] <- format(round(correl[ll], digits=max(digits-2,2)))
	  correl[!ll] <- ""
      k <- k + 1
      z[[k]] <- list(type='print',
                     list(correl[-1, - p, drop = FALSE], quote = FALSE),
                     title='Correlation Matrix for Parameter Estimates')
	}

  prModFit(x, title=title, z, digits=digits, coefs=coefs, latex=latex, ...)
}

predict.bj <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms","ccterms","adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE,
           conf.type=c('mean','individual','simultaneous'),
           kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               kint,
               na.action, expand.na, center.terms, ...)
  }

residuals.bj <- function(object, 
						 type = c("censored","censored.normalized"), ...)
{
  type <- match.arg(type)
    
  y <- object$y
  aty <- attributes(y)
  if('y' %nin% names(object)) stop('did not use y=TRUE with fit')
  ncy <- ncol(y)
  r <- y[, - ncy, drop=FALSE] - object$linear.predictors
  if(type=='censored.normalized') r <- r / object$stats['sigma']
  label(r) <- if(type=='censored') 
    'Residual' else 'Normalized Residual'
  ev <- y[, ncy]
  label(ev) <- label(y)
  units(r)  <- units(y)
  r <- Surv(r, ev)
  attr(r,'type') <- aty$type
  class(r) <- c('residuals.bj', 'Surv')
  if (length(object$na.action)) naresid(object$na.action, r)
  else r
}


validate.bj <-
  function(fit, method="boot", B=40,
           bw=FALSE, rule="aic", type="residual", sls=.05, aics=0,
           force=NULL, estimates=TRUE, pr=FALSE,
           tol=1e-7, rel.tolerance=1e-3, maxiter=15, ...)
{

  if(!(length(fit$x) && length(fit$y)))
    stop('you must specify x=TRUE and y=TRUE to bj')
  xb <- fit$linear.predictors
  ny <- dim(fit$y)
  nevents <- sum(fit$y[,ny[2]])

  ##Note: fit$y already has been transformed by the link function by psm

  distance <- function(x,y,fit,iter,evalfit=FALSE,fit.orig,
                       maxiter=15, tol=1e-7, rel.tolerance=1e-3, ...)
    {
      ##Assumes y is matrix with 1st col=time, 2nd=event indicator
      dxy.cens(x, y)["Dxy"]
    }

  predab.resample(fit, method=method,
                  fit=bj.fit2, measure=distance,
                  pr=pr, B=B, bw=bw, rule=rule, type=type,
                  sls=sls, aics=aics, force=force, estimates=estimates,
                  maxiter=maxiter, tol=tol,
                  rel.tolerance=rel.tolerance, ...)
}


bj.fit2 <- function(x,y,iter=0,maxiter=15, 
                    init=NULL, rel.tolerance=1e-3, tol=1e-7, ...)
{
  e <- y[,2]
  if(sum(e) < 1)return(list(fail=TRUE))
  x <- x	#Get around lazy evaluation creating complex expression
  f <- bj.fit(as.matrix(x), y,
              control=list(iter.max=maxiter, eps=rel.tolerance, tol=tol))
  if(f$fail) warning('bj.fit failed')
  f
}
  

latex.bj <- function(...) latexrms(...)
