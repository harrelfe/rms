orm <- function(formula, data, subset, na.action=na.delete,
				method="orm.fit", model=FALSE, x=FALSE, y=FALSE, 
				linear.predictors=TRUE, se.fit=FALSE, 
				penalty=0, penalty.matrix, tol=1e-7, eps=0.005,
                var.penalty=c('simple','sandwich'), scale=FALSE, ...)
{
  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  if(!missing(penalty) || !missing(penalty.matrix))
    stop('penalty not yet implemented')
  m <- match.call(expand.dots=FALSE)
  mc <- match(c("formula", "data", "subset", "na.action"), 
             names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  m$drop.unused.levels <- TRUE
  
  m[[1]] <- as.name("model.frame")
  nact <- NULL
  if(missing(data)) data <- NULL

  tform <- terms(formula, data=data)
  if(!missing(data) || (
						length(atl <- attr(tform,"term.labels")) && 
						any(atl!=".")))	{ ##X's present

    dul <- .Options$drop.unused.levels
    if(!length(dul) || dul) {
      on.exit(options(drop.unused.levels=dul))
      options(drop.unused.levels=FALSE)
    }

    X <- Design(eval.parent(m))
    atrx <- attributes(X)
    sformula <- atrx$sformula
    nact <- atrx$na.action
    if(method=="model.frame") return(X)
    Terms <- atrx$terms
    attr(Terms, "formula") <- formula
    atr <- atrx$Design

    Y <- model.extract(X, 'response')
    offs <- atrx$offset
    if(!length(offs)) offs <- 0
    if(model) m <- X
    X <- model.matrix(sformula, X)
    X <- X[,-1,drop=FALSE]
    dimnames(X)[[2]] <- atr$colnames
    xpres <- length(X) > 0

    p <- length(atr$colnames)
    n <- length(Y)

    penpres <- !(missing(penalty) && missing(penalty.matrix))
    if(penpres && missing(var.penalty))
      warning('default for var.penalty has changed to "simple"')

    if(!penpres) penalty.matrix <- matrix(0,ncol=p,nrow=p) else { 
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr, X) else
      if(nrow(penalty.matrix) != p || ncol(penalty.matrix) != p) stop(
             paste("penalty.matrix does not have",p,"rows and columns"))
      psetup <- Penalty.setup(atr, penalty)
      penalty <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1)
        penalty.matrix <- multiplier*penalty.matrix
      else
        {
          a <- diag(sqrt(multiplier))
          penalty.matrix <- a %*% penalty.matrix %*% a
        }
    }
  }
  else
    {
      X <- eval.parent(m)
      offs <- model.offset(X)
      if(!length(offs)) offs <- 0
      Y <- model.extract(X, 'response')
      Y <- Y[!is.na(Y)]
      Terms <- X <- NULL
      xpres <- FALSE
      penpres <- FALSE
      penalty.matrix <- NULL
    }  ##Model: y~. without data= -> no predictors
  
  if(method=="model.matrix") return(X)

  if(existsFunction(method)) {
      fitter <- getFunction(method)
      f <- fitter(X, Y, offset=offs,
                  penalty.matrix=penalty.matrix, tol=tol, eps=eps,
                  scale=scale, ...)
    }
    else stop(paste("unimplemented method:", method))
  
  if(f$fail) {
    cat("Unable to fit model using ", dQuote(method), '\n')
    return(f)
  }
  
  f$call <- NULL
  f$sformula <- sformula
  if(model) f$model <- m
  if(x) f$x <- X
  if(y) f$y <- Y
  nrp <- f$non.slopes
  if(penpres) {
    f$penalty <- penalty
    ## Get improved covariance matrix
    v <- f$var
    if(var.penalty == 'sandwich') f$var.from.info.matrix <- v
    f.nopenalty <- 
        fitter(X, Y, offset=offs, initial=f$coef, maxit=1, tol=tol)
    ##  info.matrix.unpenalized <- solvet(f.nopenalty$var, tol=tol)
    info.matrix.unpenalized <- f.nopenalty$info.matrix
    dag <- diag(info.matrix.unpenalized %*% v)
    f$effective.df.diagonal <- dag
    f$var <- if(var.penalty == 'simple') v else
       v %*% info.matrix.unpenalized %*% v
    df   <- sum(dag[-(1:nrp)])
    lr   <- f.nopenalty$stats["Model L.R."]
    pval <- 1 - pchisq(lr, df)
    f$stats[c('d.f.','Model L.R.','P')] <- c(df, lr, pval)  
  }
  ass <- if(xpres) DesignAssign(atr, nrp, Terms) else list()
  
  if(xpres) {
    if(linear.predictors) names(f$linear.predictors) <- names(Y)
    else
        f$linear.predictors <- NULL

      if(se.fit) {
        X <- cbind(1, X)
        se <- drop((((X %*% f$var) * X) %*% rep(1, ncol(X)))^.5)
        names(se) <- names(Y)
        f$se.fit <- se
      }
  }
  f <- c(f, list(call=call, Design=if(xpres)atr,
                 scale.pred=if(f$family=='logistic')
                  c("log odds","Odds Ratio") else
                  if(f$family=='loglog') c("log hazard", "Hazard Ratio"),
                 terms=Terms, assign=ass, na.action=nact))
  class(f) <- c("orm","rms")
  f
}

print.orm <- function(x, digits=4, coefs=TRUE,
                      intercepts=x$non.slopes < 10,
                      latex=FALSE,
                      title, ...) {

  if(missing(title)) {
    title <- switch(x$family,
                    logistic = 'Logistic (Proportional Odds)',
                    probit = 'Probit',
                    cauchit = 'Cauchy',
                    loglog  = '-log-log',
                    cloglog = 'Complementary log-log')
    title <- paste(title, 'Ordinal Regression Model')
  }
                    
  z <- list()
  k <- 0

  lf <- length(x$freq)
  if(lf > 3 && lf <= 20) {
    k <- k + 1
    z[[k]] <- list(type='print', list(x$freq),
                   title='Frequencies of Responses')
  }
  if(length(x$nmiss)) {  ## for backward compatibility
    k <- k + 1
    z[[k]] <- list(type='print', list(x$nmiss),
                   title='Frequencies of Missing Values Due to Each Variable')
  }
  else if(length(x$na.action)) {
    k <- k + 1
    z[[k]] <- list(type=paste('naprint', class(x$na.action), sep='.'),
                   list(x$na.action))
  }
  
  ns <- x$non.slopes
  ## coefficient intercepts kept: (fit.mult.impute)
  cik <- attr(x$coef, 'intercepts')  # esp. for fit.mult.impute
  if(length(cik) && intercepts) {
    warning('intercepts=TRUE not implemented for fit.mult.impute objects')
    intercepts <- FALSE
  }
    
  pm <- x$penalty.matrix
  penaltyFactor <- NULL
  if(length(pm)) {
    psc <- if(length(pm) == 1) sqrt(pm)
    else
      sqrt(diag(pm))
    penalty.scale <- c(rep(0, ns), psc)
    cof <- matrix(x$coef[-(1 : ns)], ncol=1)
    ## This logic does not handle fit.mult.impute objects
    k <- k + 1
    z[[k]] <- list(type='print', list(as.data.frame(x$penalty, row.names='')),
                   title='Penalty factors')
    penaltyFactor <- as.vector(t(cof) %*% pm %*% cof)
  }
  vv <- diag(vcov(x, intercepts=if(intercepts) 'all' else 'none'))
  if(!intercepts) {
    nints <- if(!length(cik)) ns else {
      if(length(cik) == 1 && cik ==0) 0 else length(cik)
    }
    ints.to.delete <- if(ns == 0 || nints == 0) integer(0) else 1:nints
    vv <- c(rep(NA, nints), vv)
  }
  cof   <- x$coef
  stats <- x$stats

  maxd <- signif(stats['Max Deriv'], 1)
  if(latex) maxd <- paste('$', latexSN(maxd), '$', sep='')
  ci <- x$clusterInfo
  misc <- reVector(Obs           = stats['Obs'],
                   'Unique Y'    = stats['Unique Y'],
                   'Cluster on'  = ci$name,
                   Clusters      = ci$n,
                   'Median Y'    = stats['Median Y'],
                   'max |deriv|' = maxd)
  if(length(x$freq) < 4) {
    names(x$freq) <- paste(if(latex)'~~' else ' ',
                           names(x$freq), sep='')
    misc <- c(misc[1], x$freq, misc[-1])
  }
  lr   <- reVector('LR chi2'    = stats['Model L.R.'],
                   'd.f.'       = round(stats['d.f.'],3),
                   'Pr(> chi2)' = stats['P'],
                   'Score chi2' = stats['Score'],
                   'Pr(> chi2)' = stats['Score P'],
                   Penalty      = penaltyFactor)
  disc <- reVector(R2=stats['R2'], g=stats['g'], gr=stats['gr'],
                   '|Pr(Y>=median)-0.5|'=stats['pdm'])
  discr <-reVector(rho=stats['rho'])
  
  headings <- list('',
                   c('Model Likelihood','Ratio Test'),
                   c('Discrimination',' Indexes'),
                   c('Rank Discrim.','Indexes'))
  data <- list(misc, c(lr, c(2,NA,-4,2,-4,2)), c(disc,3), c(discr,3))
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))

  if(coefs) {
    k <- k + 1
    if(!intercepts) {
      j   <- - ints.to.delete
      cof <- cof[j]
      vv  <- vv[j]
      if(length(pm)) penalty.scale <- penalty.scale[j]
    }
    z[[k]] <- list(type='coefmatrix',
                   list(coef=cof, se=sqrt(vv),
                        aux=if(length(pm)) penalty.scale,
                        auxname='Penalty Scale'))
  }
  
  prModFit(x, title=title, z, digits=digits,
           coefs=coefs, latex=latex, ...)
}

Mean.orm <- function(object, codes=FALSE, ...)
  Mean.lrm(object, codes=codes, ...)

Quantile.orm <- function(object, codes=FALSE, ...)
  {
    ns <- object$non.slopes
    if(ns < 5)
stop('using this function only makes sense for >5 ordered response categories')
    if(codes) vals <- 1:length(object$yunique)
    else
      {
        vals <- as.numeric(object$yunique)
        if(any(is.na(vals)))
          stop('values of response levels must be numeric for codes=FALSE')
      }
    f <- function(q=.5, lp=numeric(0), intercepts=numeric(0), values=numeric(0),
                  interceptRef=integer(0), cumprob=NULL, inverse=NULL)
      {
        ## Solve inverse(1 - q) = a + xb; inverse(1 - q) - xb = a
        ## Shift intercepts to left one position because quantile
        ## is such that Prob(Y <= y) = q whereas model is stated in
        ## Prob(Y >= y)
        lp <- lp - intercepts[interceptRef]
        ## Interpolation on linear predictors scale:
        ## z <- approx(c(intercepts, -1e100), values,
        ##             xout=inverse(1 - q) - lp, rule=2)$y
        ## Interpolation approximately on Y scale:
        lpm <- mean(lp, na.rm=TRUE)
        z <- approx(c(cumprob(intercepts + lpm), 0), values,
                    xout=cumprob(inverse(1 - q) - lp + lpm), rule=2)$y
        names(z) <- names(lp)
        z
      }
    trans <- object$trans
    formals(f) <- list(q=.5, lp=numeric(0), intercepts=object$coef[1:ns],
                       values=vals, interceptRef=object$interceptRef,
                       cumprob=trans$cumprob, inverse=trans$inverse)
    f
  }

ExProb <- function(object, ...) UseMethod("ExProb")

ExProb.orm <- function(object, codes=FALSE, ...)
  {
    ns <- object$non.slopes
    if(codes) vals <- 1:length(object$freq)
    else {
      vals <- as.numeric(object$yunique)
      if(any(is.na(vals)))
        stop('values of response levels must be numeric for codes=FALSE')
    }
    f <- function(lp=numeric(0), y=NULL, intercepts=numeric(0),
                  values=numeric(0),
                  interceptRef=integer(0), cumprob=NULL, yname=NULL)
      {
        lp <- lp - intercepts[interceptRef]
        prob <- cumprob(sapply(c(1e30, intercepts), '+', lp))
        dim(prob) <- c(length(lp), length(values))
        
        if(!length(y)) {
          colnames(prob) <- paste('Prob(Y>=', values, ')', sep='')
          return(structure(list(y=values, prob=prob, yname=yname),
                           class='ExProb'))
        }
        p <- t(apply(prob, 1,
                   function(probs) {
                     pa <- approx(values, probs, xout=y,
                                  f=1, method='constant')$y
                     pa[y < min(values)] <- 1.
                     pa[y > max(values)] <- 0.
                     pa
                   } ) )
        if(length(y) > 1) {
          colnames(p) <- paste('Prob(Y>=', y, ')', sep='')
          p <- list(y=y, prob=p)
        }
        structure(drop(p), class='ExProb')
      }
    trans <- object$trans
    formals(f) <- list(lp=numeric(0), y=NULL, intercepts=object$coef[1:ns],
                       values=vals, interceptRef=object$interceptRef,
                       cumprob=trans$cumprob, yname=all.vars(object$terms)[1])
    f
  }

plot.ExProb <- function(x, ..., data=NULL,
                        xlim=NULL, xlab=x$yname, ylab=expression(Prob(Y>=y)),
                        col=par('col'), col.vert='gray85',
                        pch=20, pch.data=21,
                        lwd=par('lwd'), lwd.data=lwd, lty.data=2, key=TRUE)
  {
    xlab <- xlab
    Y <- x[[2]]
    x <- x[[1]]
    if(!length(xlim)) xlim <- range(x)
    plot(0, 0, xlim=xlim, ylim=c(0,1), type='n', xlab=xlab, ylab=ylab, ...)
    if(!is.matrix(Y)) Y <- matrix(Y, nrow=1)
    xpts <- ypts <- numeric(0)
    for(i in 1:nrow(Y)) {
      y <- Y[i,]
      segments(x, y, x, c(y[-1], 0), col=col.vert)
      segments(x[-length(x)], y[-1], x[-1], y[-1], col=i, lwd=lwd)
      points(x, y, pch=pch, col=i)
      if(length(data)) { xpts <- c(xpts, x); ypts <- c(ypts, y) }
    }
    if(!length(data)) return(invisible())
    if(!is.list(data)) {
      groups <- rep(' ', length(data))
      Y      <- nomiss(data)
    } else {
      if(!is.list(data) || length(data) < 2) stop('inappropriate data=')
      data <- nomiss(data)
      groups <- data[[1]]
      Y      <- data[[2]]
    }
    i      <- 0
    for(g in unique(groups)) {
      i <- i + 1
      s <- groups == g
      y <- nomiss(Y[s])
      x <- sort(unique(y))
      f <- c(1., 1. - cumsum(table(y)) / length(y))
      if(x[1] > min(Y)) {
        x <- c(min(Y), x)
        f <- c(1.,     f)
      }
      y <- f[-length(f)]
      segments(x,             y,     x,     c(y[-1], 0),
               col=col.vert, lty=lty.data)
      segments(x[-length(x)], y[-1], x[-1], y[-1],
               col=i, lty=lty.data, lwd=lwd.data)
      points(x, y, pch=pch.data, col=i)
      xpts <- c(xpts, x); ypts <- c(ypts, y)
    }
    if(key && is.list(data))
      putKeyEmpty(xpts, ypts,
                  labels=unique(groups),
                  col=1:i, xlim=xlim, grid=FALSE)
    invisible()
  }
