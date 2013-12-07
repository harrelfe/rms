bootcov <- function(fit, cluster, B=200, fitter, coef.reps=TRUE, 
                    loglik=FALSE, pr=FALSE, maxit=15, eps=.0001,
                    group=NULL, stat=NULL) {

  coxcph <- inherits(fit,'coxph') || inherits(fit,'cph')

  nfit <- class(fit)[1]
  
  if(length(fit$weights) && (coxcph || nfit[1]=='Rq'))
    stop('does not handle weights')

  if(!length(X <- fit$x) | !length(Y <- fit$y))
    stop("you did not specify x=TRUE and y=TRUE in the fit")


  sc.pres <- match('scale',names(fit),0) > 0
  ns <- fit$non.slopes

  if(nfit=='psm') {
    fixed <- fit$fixed   #psm only
    fixed <- if(length(fixed) == 1 && is.logical(fixed) && !fixed) list()
    else list(scale=TRUE)
    
    fixed <-  NULL
    dist <-   fit$dist
    parms <-  fit$parms
  }

  if(nfit %in% c('Glm','orm')) fitFamily <- fit$family

  ## For orm fits, find y cutoff target for intercept (from median of
  ## original sample)
  ytarget <- if(nfit == 'orm')
    with(fit,
         ifelse(is.numeric(yunique), yunique[interceptRef + 1L],
                interceptRef + 1L))

  ## See if ordinal regression being done
  ordinal <- nfit == 'orm' || (nfit == 'lrm' && length(unique(Y)) > 2)
    
  penalty.matrix <- fit$penalty.matrix

  if(missing(fitter)) {
    fitter <-
      switch(nfit,
             ols=if(length(penalty.matrix)) {
               function(x, y, penalty.matrix,...) {
                 lm.pfit(cbind(Intercept=1., x), y,
                         penalty.matrix=penalty.matrix,
                         tol=1e-11, regcoef.only=TRUE)
               }
             }
             else function(x, y, ...) {
               lm.fit.qr.bare(x, y, tolerance=1e-11, intercept=TRUE)
             }, 
             lrm=function(x, y, maxit=15, eps=.0001, penalty.matrix,...) {
               lrm.fit(x, y, maxit=maxit, tol=1E-11, eps=eps,
                       penalty.matrix=penalty.matrix)
             }, 
             cph=function(x, y, strata=NULL, maxit=15, eps=.0001,...) {
               coxphFit(x, y, strata=strata, iter.max=maxit, 
                        eps=eps, method="efron", toler.chol=1e-11, type='right')
             },
             psm=function(x, y, maxit=15,...) {
               survreg.fit2(x, y, dist=dist, parms=parms, fixed=fixed,
                            offset=NULL, init=NULL, maxiter=maxit)
             },
             bj=function(x, y, maxit=15, eps=.0001, ...) {
               bj.fit(x, y, control=list(iter.max=maxit, eps=eps))
             },
             Glm=function(x, y, ...) {
               glm.fit(cbind(1., x), as.vector(y), family=fitFamily)
             },
             Rq=RqFit(fit, wallow=FALSE),
             orm=function(x, y, maxit=14L, eps=.005, tol=1e-7,
               ytarget=NULL, ...) {
               f <- orm.fit(x, y, family=fitFamily,
                            maxit=maxit, eps=eps, tol=tol)
               ns  <- f$non.slopes
               cof <- f$coefficients
               if(length(ytarget)) {
                 ## Y values corresponding to intercepts
                 yu <- f$yunique[-1]
                 ## Linearly interpolate to return an intercept aimed
                 ## at Y >= ytarget
                 intcept <- approx(yu, cof[1:ns], xout=ytarget)$y
                 ## if(min(abs(intcept - cof[1:ns])) > 1e-9) cat('****')
                 intattr <- approx(yu, 1:ns, xout=ytarget)$y
               }
               else {
                 k       <- f$interceptRef
                 intattr <- k
                 intcept <- cof[k]
               }
               names(intcept) <- 'Intercept'
               cof <- c(intcept, cof[(ns + 1) : length(cof)])
               attr(cof, 'intercepts') <- intattr
               f$coefficients          <- cof
               f
             })
  }
  
  if(!length(fitter)) stop("fitter not valid")
  
  if(loglik)
    {
      oosl <- switch(nfit,
                     ols=oos.loglik.ols,
                     lrm=oos.loglik.lrm,
                     cph=oos.loglik.cph,
                     psm=oos.loglik.psm,
                     Glm=oos.loglik.Glm)
      
      if(!length(oosl))
        stop('loglik=TRUE but no oos.loglik method for model in rmsMisc')
      
      Loglik <- double(B+1)
      Loglik[B+1] <- oosl(fit)
    }
  else Loglik <- NULL
  
  n     <- nrow(X)
  cof   <- fit$coefficients
  if(nfit == 'orm') {
    iref <- fit$interceptRef
    cof <- cof[c(iref, (ns + 1L) : length(cof))]
  }
  p     <- length(cof)
  vname <- names(cof)
  if(sc.pres) {
    p     <- p + 1L
    vname <- c(vname, "log scale")
  }

  ## Function to carry non-NA values backwards and replace NAs at the
  ## right end with zeros.  This will cause cell proportions for unsampled
  ## Y values to be zero for the purpose of computing means
  ## The zero placements will mess up bootstrap covariance matrix however
  fillInTheBlanks <- function(S) {
    ## http://stackoverflow.com/questions/1782704/propagating-data-within-a-vector/1783275#1783275
    ## NA in S are replaced with observed values
    ## accepts a vector possibly holding NA values and returns a vector
    ## where all observed values are carried forward and the first is
    ## also carried backward.  cfr na.locf from zoo library.
    L <- !is.na(S)
    c(S[L][1L], S[L])[cumsum(L) + 1L]
  }
  ## vn = names of full coefficient vector
  ## ns = # non-slopes (intercepts) in full vector (target)
  ## nc = # non-slopes for current fit in cof
  fill <- function(cof, vn, ns) {
    p <- length(vn)
    if(length(cof) == p) return(cof)
    nc <- ns - (p - length(cof))
    cints              <- cof[1L : nc]  ## current intercepts
    ints               <- rep(NA, ns)
    names(ints)        <- vn[1L : ns]
    ints[names(cints)] <- cints
    ## Set not just last intercept to -Inf if missing but set all
    ## NA intercepts at the right end to -Inf.  This will later lead to
    ## cell probabilities of zero for bootstrap-omitted levels of Y
    if(is.na(ints[ns])) {
      l <- ns
      if(ns > 1L) {
        for(j in (ns - 1L) : 1L) {
          if(!is.na(ints[j])) break
          l <- j
        }
      }
      ints[l : ns] <- -Inf  ## probability zero of exceeding unobserved high Y
    }
    #### CHANGE TO FILL IN ONLY INTERCEPTS
    c(rev(fillInTheBlanks(rev(ints))), cof[-(1L : nc)])
  }
    
  bar <- rep(0, p)
  cov <- matrix(0, nrow=p, ncol=p, dimnames=list(vname,vname))
  if(coef.reps) coefs <- matrix(NA, nrow=B, ncol=p, dimnames=list(NULL,vname))
  if(length(stat)) stats <- numeric(B)
  
  Y <- as.matrix(if(is.factor(Y)) unclass(Y) else Y)
  ny <- ncol(Y)

  Strata <- fit$Strata
  
  nac <- fit$na.action
  
  if(length(group)) {
    if(length(group) > n) {
      ## Missing observations were deleted during fit
      if(length(nac)) {
        j <- !is.na(naresid(nac, Y) %*% rep(1,ny))
        group <- group[j]
      }
    }
    
    if(length(group) != n)
      stop('length of group does not match # rows used in fit')

    group.inds <- split(1:n, group)
    ngroup     <- length(group.inds)
  }
  else ngroup <- 0

  anyinf <- FALSE

  if(!exists('.Random.seed')) runif(1)
  seed <- .Random.seed
  if(missing(cluster)) {
    b <- 0
    pb <- setPb(B, type='Boot', onlytk=!pr, every=20)
    for(i in 1:B) {
      pb(i)
      
      if(ngroup) {
        j <- integer(n)
        for(si in 1L : ngroup) {
          gi    <- group.inds[[si]]
          j[gi] <- sample(gi, length(gi), replace=TRUE)
        }
      }
      else j <- sample(1L : n, n, replace=TRUE)
      
      ## Note: If Strata is NULL, NULL[j] is still NULL
      
      f <- tryCatch(fitter(X[j,,drop=FALSE], Y[j,,drop=FALSE], maxit=maxit, 
                           eps=eps, ytarget=ytarget,
                           penalty.matrix=penalty.matrix, strata=Strata[j]),
                    error=function(...) list(fail=TRUE))
      if(length(f$fail) && f$fail) next
          
      cof <- f$coefficients
      if(any(is.na(cof))) next   # glm
      b <- b + 1L
          
      if(sc.pres) cof <- c(cof, 'log scale' = log(f$scale))

      ## Index by names used since some intercepts may be missing in a
      ## bootstrap resample from an ordinal logistic model
      ## Missing coefficients represent values of Y not appearing in the
      ## bootstrap sample.  Carry backwards the next non-NA intercept
      if(ordinal) cof <- fill(cof, vname, ns)
      if(any(is.infinite(cof))) anyinf <- TRUE
      if(coef.reps) coefs[b,] <- cof
      
      if(length(stat)) stats[b] <- f$stats[stat]
      
      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)
      
      if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
    }
  }
  else {
    if(length(cluster) > n) {
      ## Missing obs were deleted during fit
      if(length(nac)) {
        j <- !is.na(naresid(nac, Y) %*% rep(1,ny))
        cluster <- cluster[j]
      }
    }
    
    if(length(cluster)!=n)
      stop("length of cluster does not match # rows used in fit")
    
    if(any(is.na(cluster))) stop("cluster contains NAs")
    
    cluster <- as.character(cluster)
    
    clusters <- unique(cluster)
    nc <- length(clusters)
    Obsno <- split(1:n, cluster)
    
    b <- 0
    pb <- setPb(B, type='Boot', onlytk=!pr, every=20)

    for(i in 1L : B) {
      pb(i)
      
      ## Begin addition Bill Pikounis
      if(ngroup) {
        j <- integer(0L)
        for(si in 1L : ngroup) {
          gi <- group.inds[[si]]
          cluster.gi <- cluster[gi]
          clusters.gi <- unique(cluster.gi)
          nc.gi <- length(clusters.gi)
          Obsno.gci <- split(gi, cluster.gi)
          j.gci <- sample(clusters.gi, nc.gi, replace = TRUE)
          obs.gci <- unlist(Obsno.gci[j.gci])
          j <- c(j, obs.gci)
        }
        obs <- j
      }
      else {
        ## End addition Bill Pikounis (except for closing brace below)
        j <- sample(clusters, nc, replace=TRUE)
        obs <- unlist(Obsno[j])
      }
      
      f <- tryCatch(fitter(X[obs,,drop=FALSE], Y[obs,,drop=FALSE], 
                           maxit=maxit, eps=eps, ytarget=ytarget,
                           penalty.matrix=penalty.matrix,
                           strata=Strata[obs]),
                    error=function(...) list(fail=TRUE))
      if(length(f$fail) && f$fail) next
      
      cof <- f$coefficients
      if(any(is.na(cof))) next  # glm
      b <- b + 1L
          
      if(sc.pres) cof <- c(cof, 'log scale' = log(f$scale))
      
      cof <- fill(cof, vname, ns)
      if(any(is.infinite(cof))) anyinf <- TRUE
      if(coef.reps)    coefs[b,] <- cof
      if(length(stat)) stats[b] <- f$stats[stat]
      
      bar <- bar + cof
      cof <- as.matrix(cof)
      cov <- cov + cof %*% t(cof)
      if(loglik) Loglik[b] <- oosl(f, matxv(X,cof), Y)
    }
  }
  
  if(b < B) {
    warning(paste('fit failure in',B-b,
                  'resamples.  Might try increasing maxit'))
    if(coef.reps) coefs <- coefs[1L : b,,drop=FALSE]
    Loglik <- Loglik[1L : b]
  }
  if(nfit == 'orm') attr(coefs, 'intercepts') <- iref
  
  if(anyinf) warning('at least one resample excluded highest Y values, invalidating bootstrap covariance matrix estimate')
  
  bar           <- bar / b
  fit$B         <- b
  fit$seed      <- seed
  names(bar)    <- vname
  fit$boot.coef <- bar
  if(coef.reps) fit$boot.Coef <- coefs
  
  bar <- as.matrix(bar)
  cov <- (cov - b * bar %*% t(bar)) / (b - 1L)
  fit$orig.var <- fit$var
  fit$var <- cov
  fit$boot.loglik <- Loglik
  if(length(stat)) fit$boot.stats <- stats
  if(nfit=='Rq') {
    newse <- sqrt(diag(cov))
    newt <- fit$summary[, 1L]/newse
    newp <- 2. * (1. - pt(abs(newt), fit$stats['n'] - fit$stats['p']))
    fit$summary[, 2L : 4L] <- cbind(newse, newt, newp)
  }
  fit
}
  
bootplot <- function(obj, which, X,
                     conf.int=c(.9,.95,.99),
                     what=c('density','qqnorm'),
                     fun=function(x)x,
                     labels., ...) {

  what <- match.arg(what)
  Coef <- obj$boot.Coef
  if(length(Coef)==0) stop('did not specify "coef.reps=TRUE" to bootcov')
  
  if(missing(which)) {
    if(!is.matrix(X)) X <- matrix(X, nrow=1)
    
    qoi <- matxv(X, Coef, bmat=TRUE)  # X %*% t(Coef)   ##nxp pxB = nxB
    if(missing(labels.)) {
      labels. <- dimnames(X)[[1]]
      if(length(labels.)==0) {
        labels. <- as.character(1:nrow(X))
      }
    }
  }
  else {
    qoi <- t(Coef[, which, drop=FALSE])
    nns <- num.intercepts(obj)
    if(missing(labels.)) {
      labels. <- paste(ifelse(which > nns, 'Coefficient of ', ''), 
                       dimnames(Coef)[[2]][which], sep='')
    }
  }
  
  nq <- nrow(qoi)
  qoi <- fun(qoi)
  quan <- NULL
  
  if(what=='density') {
    probs <- (1+conf.int)/2
    probs <- c(1-probs, probs)
    quan <- matrix(NA, nrow=nq, ncol=2*length(conf.int),
                   dimnames=list(labels., format(probs)))

      for(j in 1:nq) {
        histdensity(qoi[j,], xlab=labels.[j], ...)
        quan[j,] <- quantile(qoi[j,], probs, na.rm=TRUE)
        abline(v=quan[j,], lty=2)
        title(sub=paste('Fraction of effects>',fun(0),' = ',
                format(mean(qoi[j,]>fun(0))),sep=''),adj=0)
      }
  }
  else {
    for(j in 1:nq) {
      qqnorm(qoi[j,], ylab=labels.[j])
      qqline(qoi[j,])
    }
  }
  
  invisible(list(qoi=drop(qoi), quantiles=drop(quan)))
}


## histdensity runs hist() and density(), using twice the number of
## class than the default for hist, and 1.5 times the width than the default
## for density

histdensity <- function(y, xlab, nclass, width, mult.width=1, ...) {
  y <- y[is.finite(y)]
  if(missing(xlab)) {
    xlab <- label(y)
    if(xlab=='') xlab <- as.character(sys.call())[-1]
  }

  if(missing(nclass)) nclass <- (logb(length(y),base=2)+1)*2
  
  hist(y, nclass=nclass, xlab=xlab, probability=TRUE, ...)
  if(missing(width)) {
    nbar <- logb(length(y), base = 2) + 1
    width <- diff(range(y))/nbar*.75*mult.width
  }

  lines(density(y,width=width,n=200))
  invisible()
}


confplot <- function(obj, X, against, 
                     method=c('simultaneous', 'pointwise'),
                     conf.int=0.95,
                     fun=function(x) x, 
                     add=FALSE, lty.conf=2, ...) {

  method <- match.arg(method)
  if(length(conf.int)>1) stop('may not specify more than one conf.int value')

  boot.Coef <- obj$boot.Coef
  if(length(boot.Coef)==0) stop('did not specify "coef.reps=TRUE" to bootcov')
  
  if(!is.matrix(X)) X <- matrix(X, nrow=1)
  
  fitted <- fun(matxv(X, obj$coefficients))
  
  if(method=='pointwise') {
    pred <- matxv(X, boot.Coef, bmat=TRUE)   ## n x B
    p <- fun(apply(pred, 1, quantile,
                   probs=c((1 - conf.int)/2, 1 - (1 - conf.int)/2),
                   na.rm=TRUE))
    lower <- p[1,]
    upper <- p[2,]
  }
  else {
    boot.Coef <- rbind(boot.Coef, obj$coefficients)
    loglik    <- obj$boot.loglik
    if(length(loglik)==0) stop('did not specify "loglik=TRUE" to bootcov')
    
    crit  <- quantile(loglik, conf.int, na.rm=TRUE)
    qual  <- loglik <= crit
    boot.Coef <- boot.Coef[qual,,drop=FALSE]
    pred   <- matxv(X, boot.Coef, bmat=TRUE)  ## n x B
    upper  <- fun(apply(pred, 1, max))
    lower  <- fun(apply(pred, 1, min))
    pred   <- fun(pred)
  }
  
  if(!missing(against)) {
    lab <- label(against)
    if(lab=='') lab <- (as.character(sys.call())[-1])[3]
    
    if(add) lines(against, fitted, ...)
    else plot(against, fitted, xlab=lab, type='l', ...)
    
    lines(against, lower, lty=lty.conf)
    lines(against, upper, lty=lty.conf)
  }
  if(missing(against)) list(fitted=fitted, upper=upper, lower=lower)
  else invisible(list(fitted=fitted, upper=upper, lower=lower))
}

# Construct object suitable for boot:boot.ci
# Use boot package to get BCa confidence limits for a linear combination of
# model coefficients, e.g. bootcov results boot.Coef
# If boot.ci fails return only ordinary percentile CLs
bootBCa <- function(estimate, estimates, type=c('percentile','bca','basic'),
                    n, seed, conf.int=0.95) {
  type <- match.arg(type)
  if(type != 'percentile' && !require(boot)) stop('boot package not installed')
  estimate <- as.vector(estimate)
  ne <- length(estimate)
  if(!is.matrix(estimates)) estimates <- as.matrix(estimates)
  if(ncol(estimates) != ne)
    stop('no. columns in estimates != length of estimate')
  if(type == 'percentile') {
    a <- apply(estimates, 2, quantile,
                 probs=c((1-conf.int)/2, 1-(1-conf.int)/2), na.rm=TRUE)
    if(ne == 1) a <- as.vector(a)
    return(a)
  }
  lim <- matrix(NA, nrow=2, ncol=ne, dimnames=list(c('Lower','Upper'),NULL))
  R <- nrow(estimates)
  for(i in 1:ne) {
    w <- list(sim= 'ordinary',
              stype = 'i',
              t0 = estimate[i],
              t  = estimates[,i,drop=FALSE],
              R  = R,
              data    = 1:n,
              strata  = rep(1,   n),
              weights = rep(1/n, n),
              seed = seed,
              statistic = function(...) 1e10,
              call = match.call())

    cl <- try(boot.ci(w, type=type, conf=conf.int), silent=TRUE)
    if(inherits(cl, 'try-error')) {
      cl <- c(NA,NA)
    warning('could not obtain bootstrap confidence interval')
    } else {
      cl <- if(type == 'bca') cl$bca else cl$basic
      m <- length(cl)
      cl <- cl[c(m - 1, m)]
    }
    lim[,i] <- cl
  }
  if(ne==1) as.vector(lim) else lim
}
