pentrace <-
  function(fit, penalty, penalty.matrix,
					 method=c('grid', 'optimize'),
					 which=c('aic.c', 'aic', 'bic'), target.df=NULL,
           fitter, pr=FALSE,
           tol=1e-13, keep.coef=FALSE, complex.more=TRUE,
           verbose=FALSE, maxit=20, subset, noaddzero=FALSE, ...)
{
  ## Need to check Strata for cph

  method <- match.arg(method)
  which  <- match.arg(which)
  tdf    <- length(target.df)
  if(tdf) method <- 'optimize'

  if(! length(X <- fit$x) || ! length(Y <- as.matrix(fit$y)))
    stop("you did not specify x=TRUE and y=TRUE in the fit")
  fit$x <- fit$y <- NULL

##  if(length(pn <- fit$penalty) > 0 && max(unlist(pn)) != 0)
##    warning('you did not specify penalty= in fit so that unpenalized model can be a candidate for the best model')

  sc.pres <- match("parms", names(fit), 0) > 0

  fixed <- NULL
  dist  <- fit$dist
  parms <- fit$parms

  clas <- class(fit)[1]
  isols <- clas=='ols'

  if(!(isols || inherits(fit, 'lrm')))
    stop("at present pentrace only works for lrm or ols")

  if(missing(fitter))
    fitter <- switch(clas,
                     ols=function(x, y, maxit, ...) lm.pfit(x, y, ...),
                     lrm=function(x, y, maxit=20, ...)
                            lrm.fit(x, y, maxit=maxit, compvar = TRUE,  ...),
                     cph=function(x, y, maxit=20, ...) coxphFit(x, y,
                                        strata=Strata, iter.max=maxit,
                                        eps=.0001, method="efron",
                                        toler.chol=tol),
                     psm=function(x, y, maxit=20,...) survreg.fit2(x, y,
                                        dist=dist, parms=parms, fixed=fixed,
                                        offset=NULL,
                                        init=NULL, maxiter=maxit))

  if(!length(fitter))stop("fitter not valid")

  Strata <- fit$strata

  if(!missing(subset)) {
    Y <- Y[subset,, drop=FALSE]
    X <- X[subset,, drop=FALSE]
    Strata <- Strata[subset, drop=FALSE]  # NULL[] is NULL
  }
  n <- nrow(Y)

  atr <- fit$Design

  if(missing(penalty.matrix))
    penalty.matrix <- Penalty.matrix(atr, X)

  obj.best <- -1e10

  ns <- num.intercepts(fit)


  islist <- is.list(penalty)

  if(islist) {
    penalty <- expand.grid(penalty)
    if(complex.more && ncol(penalty) > 1 && nrow(penalty) > 1) {
      ikeep <- NULL
      for(i in 1:nrow(penalty)) {
        ok <- TRUE
        peni <- penalty[i,]
        for(j in 2:length(peni))
          if(peni[[j]] < peni[[j-1]]) ok <- FALSE
        if(ok) ikeep <- c(ikeep, i)
      }
      penalty <- penalty[ikeep,,drop=FALSE]
    }
    np <- nrow(penalty)
  }
  else {
    if(method == 'grid' && ! noaddzero) penalty <- c(0, penalty[penalty > 0])
    np <- length(penalty)
  }

  if(method=='optimize') {
    stop('method="optimize" not yet implemented in R')

    if((islist && nrow(penalty) > 1) || (!islist && length(penalty) > 1))
      stop('may not specify multiple potential penalties when method="optimize"')

    objective <- function(pen, X, Y, z) {

      ##Problem with sending so many auxiliary parameters to nlminb -
      ##nlminb's internal parameters got shifted somehow
      n <- z$n; penalty.matrix <- z$penalty.matrix; pennames <- z$pennames
      isols <- z$isols; islist <- z$islist; tol <- z$tol; maxit <- z$maxit
      ns <- z$ns; fitter <- z$fitter; pr <- z$pr; atr <- z$atr;
      tdf <- length(z$target.df)

      if(length(pen) > 1) {
        pen <- structure(as.list(pen), names=pennames)
        penfact <- Penalty.setup(atr, pen)$multiplier
      } else penfact <- pen

      if(length(penfact)==1 || !islist) pm <- penfact*penalty.matrix
      else {
        a <- diag(sqrt(penfact))
        pm <- a %*% penalty.matrix %*% a
      }
      f <- fitter(X, Y, penalty.matrix=pm, tol=tol, maxit=maxit, ...)
      if(length(f$fail) && f$fail)
        stop('fitter failed.  Try changing maxit or tol')

      if(isols) {
        ## ols (from lm.pfit) already stored correct LR chisq and effective df
        stats <- f$stats
        df <- stats['d.f.']
        lr <- stats['Model L.R.']
        dag <- f$effective.df.diagonal
      }
      else  {
        v <- f$var   #Later: vcov(f)
        f.nopenalty <- fitter(X, Y, initial=f$coef, maxit=1, tol=tol, ...)
        if(length(f.nopenalty$fail) && f.nopenalty$fail)
          stop('fitter failed.  Try changing tol')
        info.matrix.unpenalized <-
          if(length(f.nopenalty$info.matrix))
            f.nopenalty$info.matrix else
        solvet(f.nopenalty$var, tol=tol) # -> vcov
        dag <- diag(info.matrix.unpenalized %*% v)
        df <- if(ns==0)sum(dag) else sum(dag[-(1:ns)])
        lr <- f.nopenalty$stats["Model L.R."]
      }
      obj <- switch(z$which,
                    aic.c <- lr - 2*df*(1 + (df + 1) / (n - df - 1)),
                    aic   <- lr - 2 * df,
                    bic   <- lr - df * logb(n))
      if(tdf) obj <- abs(df - z$target.df)
      if(pr) {
        w <- if(tdf) df else obj
        names(w) <- NULL
        pp <- if(islist) unlist(pen) else c(Penalty=pen)
        print(c(pp, Objective=w))
      }
      if(!tdf) obj <- -obj
      else
        attr(obj,'df') <- df
      obj
    }
    res <- nlminb(unlist(penalty), objective, lower=0, X=X, Y=Y,
                  z=list(n=n,
                    penalty.matrix=penalty.matrix, pennames=names(penalty),
                    isols=isols, islist=islist, tol=tol, maxit=maxit, ns=ns,
                    fitter=fitter, atr=atr, pr=pr, which=which,
                    target.df=target.df),
                  control=list(abs.tol=.00001,
                    rel.tol=if(tdf)1e-6 else .00001))
    return(list(penalty=if(islist)
                structure(as.list(res$parameters),names=names(penalty))
    else res$parameters,
                objective=if(tdf)res$aux$df else -res$objective))
  }

  df <- aic <- bic <- aic.c <-
    if(islist) double(length(penalty[[1]])) else double(length(penalty))

  for(i in 1 : np) {
    if(islist) {
      pen     <- penalty[i,]
      penfact <- Penalty.setup(atr, pen)$multiplier
    } else {
      pen     <- penalty[i]
      penfact <- pen
    }
    unpenalized <- all(penfact==0)

    if(i==1) Coef <- if(keep.coef) matrix(NA,ncol=length(fit$coef),nrow=np)
    else NULL

   if(unpenalized) f <- fit
    else {
      if(length(penfact) == 1 || !islist) pm <- penfact * penalty.matrix
      else {
        a <- diag(sqrt(penfact))
        pm <- a %*% penalty.matrix %*% a
      }
      f <- fitter(X, Y, penalty.matrix=pm, tol=tol, maxit=maxit, ...)
      if(length(f$fail) && f$fail)
        stop('fitter failed.  Try changing maxit or tol')
    }

    if(keep.coef) Coef[i,] <- f$coef

    if(unpenalized || isols) {
      ## ols (from lm.pfit) already stored correct LR chisq and effective df
      stats <- f$stats
      df[i] <- stats['d.f.']
      lr    <- stats['Model L.R.']
      dag <- if(unpenalized) rep(1, length(df[i]))
      else
        f$effective.df.diagonal
    }
    else {
      v <- f$var   #Later: vcov(f, regcoef.only=T)
      f.nopenalty <- fitter(X, Y, initial=f$coef, maxit=1, tol=tol, ...)
      if(length(f.nopenalty$fail) && f.nopenalty$fail)
        stop('fitter failed.  Try changing tol')
      info.matrix.unpenalized <-
        if(length(f.nopenalty$info.matrix)) f.nopenalty$info.matrix
        else
          solvet(f.nopenalty$var, tol=tol) # -> vcov
      dag <- diag(info.matrix.unpenalized %*% v)
      df[i] <- if(ns == 0)sum(dag) else sum(dag[- (1 : ns)])
      lr <- f.nopenalty$stats["Model L.R."]
      if(verbose) {
        cat('non slopes',ns,'\neffective.df.diagonal:\n')
        print(dag)
      }
    }
    aic[i]   <- lr - 2 * df[i]
    bic[i]   <- lr - df[i] * logb(n)
    aic.c[i] <- lr - 2 * df[i] * (1 + (df[i] + 1) / (n - df[i] - 1))
    obj <- switch(which, aic.c=aic.c[i], aic=aic[i], bic=bic[i])

    if(obj > obj.best) {
      pen.best <- pen
      df.best <- df[i]
      obj.best <- obj
      f.best <- f
      var.adj.best <- if(unpenalized || isols) f$var
      else
        v %*% info.matrix.unpenalized %*% v
      diag.best <- dag
    }
    if(pr) {
      d <- if(islist) as.data.frame(pen, row.names='') else
      data.frame(penalty=pen, row.names='')
      d$df <- df[i]
      d$aic <- aic[i]
      d$bic <- bic[i]
      d$aic.c <- aic.c[i]
      cat('\n'); print(d)
    }
  }
  mat <- if(islist) as.data.frame(penalty)
  else
    data.frame(penalty=penalty)
  mat$df    <- df
  mat$aic   <- aic
  mat$bic   <- bic
  mat$aic.c <- aic.c

  structure(list(penalty=pen.best, df=df.best, objective=obj.best,
                 fit=f.best, var.adj=var.adj.best, diag=diag.best,
                 results.all=mat, Coefficients=Coef), class="pentrace")
}

plot.pentrace <- function(x, method=c('points', 'image'),
						  which=c('effective.df', 'aic', 'aic.c', 'bic'),
						  pch=2, add=FALSE, ylim, ...)
{
  method <- match.arg(method)

  x     <- x$results.all

  penalty      <- x[[1]]
  effective.df <- x$df
  aic          <- x$aic
  bic          <- x$bic
  aic.c        <- x$aic.c

  if(length(x) == 5) {  ## only one variable given to penalty=

    if('effective.df' %in% which) {
      if(add) lines(penalty, effective.df) else
      plot(penalty, effective.df, xlab="Penalty", ylab="Effective d.f.",
           type="l", ...)
      if(length(which) == 1) return(invisible())
    }

    if(!add) plot(penalty, aic,
                  ylim=if(missing(ylim)) range(c(aic, bic)) else ylim,
                  xlab="Penalty",
                  ylab=expression(paste("Information Criterion (", chi^2,
                      " scale)")),
                  type=if('aic' %in% which)"l" else "n", lty=3, ...)
    else
      if('aic' %in% which) lines(penalty, aic,   lty=3, ...)
    if('bic'   %in% which) lines(penalty, bic,   lty=2, ...)
    if('aic.c' %in% which) lines(penalty, aic.c, lty=1, ...)
    if(!add && length(setdiff(which, 'effective.df')) > 1)
      title(sub=paste(if('aic.c' %in% which) "Solid: AIC_c",
              if('aic' %in% which) "Dotted: AIC",
              if('bic' %in% which) "Dashed: BIC",sep='  '),
            adj=0,cex=.75)

    return(invisible())
  }

  ## At least two penalty factors
  if(add) stop('add=TRUE not implemented for >=2 penalty factors')

  X1 <- x[[1]]
  X2 <- x[[2]]
  nam <- names(x)
  x1 <- sort(unique(X1))
  x2 <- sort(unique(X2))
  n1 <- length(x1)
  n2 <- length(x2)

  aic.r <- rank(aic); aic.r <- aic.r/max(aic.r)

  if(method=='points') {
    plot(0, 0, xlim=c(1,n1), ylim=c(1,n2), xlab=nam[1], ylab=nam[2],
         type='n', axes=FALSE, ...)
    mgp.axis(1, at=1:n1, labels=format(x1))
    mgp.axis(2, at=1:n2, labels=format(x2))
    ix <- match(X1, x1)
    iy <- match(X2, x2)
    for(i in 1:length(aic)) points(ix[i], iy[i], pch=pch, cex=(.1+aic.r[i])*3)
    return(invisible(aic.r))
  }

  z <- matrix(NA,nrow=n1,ncol=n2)
  for(i in 1:n1)
    for(j in 1:n2) z[i,j] <- aic.r[X1==x1[i] & X2==x2[j]]
  image(x1, x2, z, xlab=nam[1], ylab=nam[2], zlim=range(aic.r), ...)
  invisible(aic.r)
}

print.pentrace <- function(x, ...)
{
  cat('\nBest penalty:\n\n')
  pen <- if(is.list(x$penalty)) as.data.frame(x$penalty,row.names='')
  else
    data.frame(penalty=x$penalty, row.names='')
  pen$df <- x$df
  pen$aic <- x$aic
  print(pen)
  cat('\n')
  if(is.data.frame(x$results.all)) print(x$results.all, row.names=FALSE) else
  print(as.data.frame(x$results.all,), row.names=FALSE)
#                      row.names=rep('',length(x$results.all[[1]]))))
  invisible()
}

effective.df <- function(fit, object)
{
  atr <- fit$Design

  dag <- if(missing(object)) fit$effective.df.diagonal
  else
    object$diag
  if(length(dag)==0) stop('object not given or fit was not penalized')

  ia.or.nonlin <- param.order(atr, 2)
  nonlin       <- param.order(atr, 3)
  ia           <- param.order(atr, 4)
  ia.nonlin    <- param.order(atr, 5)

  ns <- num.intercepts(fit)
  if(ns > 0) dag <- dag[-(1:ns)]

  z <- rbind(c(length(dag),        sum(dag)),
             c(sum(!ia.or.nonlin), sum(dag[!ia.or.nonlin])),
             c(sum(ia.or.nonlin),  sum(dag[ia.or.nonlin])),
             c(sum(nonlin),        sum(dag[nonlin])),
             c(sum(ia),            sum(dag[ia])),
             c(sum(ia.nonlin),     sum(dag[ia.nonlin])))

  dimnames(z) <- list(c('All','Simple Terms','Interaction or Nonlinear',
                        'Nonlinear', 'Interaction','Nonlinear Interaction'),
                      c('Original','Penalized'))

  cat('\nOriginal and Effective Degrees of Freedom\n\n')
  print(round(z,2))
  invisible(z)
}
