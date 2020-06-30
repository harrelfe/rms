contrast <- function(fit, ...) UseMethod("contrast")

contrast.rms <-
  function(fit, a, b, a2, b2, cnames=NULL, fun=NULL, funint=TRUE,
           type=c('individual','average','joint'),
           conf.type=c('individual','simultaneous'), usebootcoef=TRUE,
           boot.type=c('percentile','bca','basic'),
           posterior.summary=c('mean', 'median', 'mode'),
           weights='equal', conf.int=0.95, tol=1e-7, expand=TRUE, ...)
{
  type <- match.arg(type)
  conf.type <- match.arg(conf.type)
  boot.type <- match.arg(boot.type)
  posterior.summary <- match.arg(posterior.summary)

  draws <- fit$draws
  bayes <- length(draws) > 0

  if(bayes & (type == 'joint' || conf.type == 'simultaneous'))
    stop('type=joint or conf.type=simultaneous not allowed for Bayesian models')
  
  zcrit <- if(length(idf <- fit$df.residual)) qt((1 + conf.int) / 2, idf) else
              qnorm((1 + conf.int) / 2)
  bcoef <- if(usebootcoef) fit$boot.Coef

  pmode <- function(x) {
    dens <- density(x)
    dens$x[which.max(dens$y)[1]]
    }
  
  if(! bayes) betas <- coef(fit)
  fite  <- fit
  if(inherits(fit, 'orm')) {
    nrp <- 1
    ## Note: is 1 for orm because vcov defaults to intercepts='mid'
    w <- c(fit$interceptRef, (num.intercepts(fit) + 1) : length(betas))
    betas <- betas[w]
    fite$coefficients <- betas    # for simult confint
    if(usebootcoef) bcoef <- bcoef[, w, drop=FALSE]
  } else nrp <- num.intercepts(fit, 'var')
    
  if(length(bcoef) && conf.type != 'simultaneous')
    conf.type <- switch(boot.type,
                        percentile = 'bootstrap nonparametric percentile',
                        bca        = 'bootstrap BCa',
                        basic      = 'basic bootstrap')
  
  da <- do.call('gendata', list(fit, factors=a, expand=expand))
  xa <- predict(fit, da, type='x')
  ma <- nrow(xa)

  if(missing(b)) {
    xb <- 0*xa
    db <- da
  } else {
    db <- do.call('gendata', list(fit, factors=b, expand=expand))
    xb <- predict(fit, db, type='x')
  }
  mb <- nrow(xb)

  if(! missing(a2)) {
    if(missing(b) || missing(b2)) stop('b and b2 must be given if a2 is given')
    da2 <- do.call('gendata', list(fit, factors=a2, expand=expand))
    xa2 <- predict(fit, da2, type='x')
    ma2 <- nrow(xa2)
    db2 <- do.call('gendata', list(fit, factors=b2, expand=expand))
    xb2 <- predict(fit, db2, type='x')
    mb2 <- nrow(xb2)
  }

  allsame <- function(x) diff(range(x)) == 0
  
  vary <- NULL
  mall <- c(ma, mb)
  ncols <- c(ncol(da), ncol(db))
  if(! missing(a2)) {
    mall <- c(mall, ma2, mb2)
    ncols <- c(ncols, ncol(da2), ncol(db2))
  }
  
  if(allsame(mall) && ! allsame(ncols)) stop('program logic error')
  if(any(sort(names(da)) != sort(names(db))))
    stop('program logic error')
  if(! missing(a2) && (any(sort(names(da)) != sort(names(da2))) ||
                       any(sort(names(da)) != sort(names(db2)))))
    stop('program logic error')
    
  if(type != 'average' && ! length(cnames)) {
    ## If all lists have same length, label contrasts by any variable
    ## that has the same length and values in all lists
    k <- integer(0)
    nam <- names(da)
    for(j in 1 : length(da)) {
      w <- nam[j]
      eq <- all(as.character(da[[w]]) == as.character(db[[w]]))
      if(! missing(a2))
        eq <- eq & all(as.character(da[[w]]) == as.character(da2[[w]])) &
          all(as.character(da[[2]]) == as.character(db2[[w]]))
      if(eq) k <- c(k, j)
    }
    if(length(k)) vary <- da[k]
  } else if(max(mall) > 1) {
    ## Label contrasts by values of longest variable in list if
    ## it has the same length as the expanded design matrix
    d <- if(ma > 1) a else b
    if(! missing(a2) && (max(ma2, mb2) > max(ma, mb)))
      d <- if(ma2 > 1) a2 else b2
    l <- sapply(d, length)
    vary <- if(sum(l == max(mall)) == 1) d[l == max(mall)]
  }

  if(sum(mall > 1) > 1 && ! allsame(mall[mall > 1]))
    stop('lists of settings with more than one row must all have the same # rows')
  mm <- max(mall)
  if(mm > 1 && any(mall == 1)) {
    if(ma == 1) xa <- matrix(xa, nrow=mm, ncol=ncol(xa), byrow=TRUE)
    if(mb == 1) xb <- matrix(xb, nrow=mm, ncol=ncol(xb), byrow=TRUE)
    if(! missing(a2)) {
      if(ma2 == 1) xa2 <- matrix(xa2, nrow=mm, ncol=ncol(xa2), byrow=TRUE)
      if(mb2 == 1) xb2 <- matrix(xb2, nrow=mm, ncol=ncol(xb2), byrow=TRUE)
    }
  }

  if(bayes && length(fun) && inherits(fit, 'blrm')) {
    if(! missing(a2)) stop('fun= is only implemented for blrm fits')
    if(missing(b))    stop('b must be specified when fun= is given')

    nd   <- nrow(draws)
    nrp  <- num.intercepts(fit)
    kint <- fit$interceptRef
    pa   <- t(matxv(xa, draws, kint=kint,  bmat=TRUE))
    pb   <- t(matxv(xb, draws, kint=kint,  bmat=TRUE))
    if(length(cnames)) colnames(pa) <- colnames(pb) <- cnames
    # If fun has an intercepts argument, the intecept vector must be
    # updated for each draw
    for(i in 1 : nd) {
      pa[i, ] <- if(funint)
                   fun(pa[i, ], intercepts=draws[i, 1L : nrp])
                 else
                   fun(pa[i, ])
      pb[i, ] <- if(funint)
                   fun(pb[i, ], intercepts=draws[i, 1L : nrp])
                 else
                   fun(pb[i, ])
    }
    if(! length(cnames))
      cnames <- if(length(vary)) rep('', ncol(pa)) else
                          as.character(1 : ncol(pa))
    colnames(pa) <- colnames(pb) <- cnames

    res <- list(esta=pa, estb=pb,
                Xa=xa, Xb=xb, 
                nvary=length(vary))
    return(structure(res, class='contrast.rms'))
    }
  
  X <- xa - xb
  if(! missing(a2)) X <- X - (xa2 - xb2)
  m <- nrow(X)
  if(nrp > 0) X <- cbind(matrix(0., nrow=m, ncol=nrp), X)
  
  if(is.character(weights)) {
    if(weights != 'equal') stop('weights must be "equal" or a numeric vector')
    weights <- rep(1,  m)
  } else if(length(weights) > 1 && type != 'average')
      stop('can specify more than one weight only for type="average"')
    else if(length(weights) != m) stop(paste('there must be', m, 'weights'))
  weights <- as.vector(weights)
  if(m > 1 && type=='average')
    X <- matrix(apply(weights*X, 2, sum) / sum(weights), nrow=1,
                dimnames=list(NULL, dimnames(X)[[2]]))

  cdraws <- NULL
  
  if(bayes) {
    cdraws <- draws %*% t(X)
    if(length(cnames)) colnames(cdraws) <- cnames
    v      <- var(cdraws)
    ndf    <- if(is.matrix(v)) nrow(v) else 1
    ci     <- apply(cdraws, 2, HPDint, prob=conf.int)
    lower  <- ci[1, ]
    upper  <- ci[2, ]
    PP     <- apply(cdraws, 2, function(u) mean(u > 0))
    se     <- apply(cdraws, 2, sd)

    est    <- switch(posterior.summary,
                     mode   = apply(cdraws, 2, pmode),
                     mean   = colMeans(cdraws),
                     median = apply(cdraws, 2, median))
    P <- Z <- NULL
  }
  else {
  est <- matxv(X, betas)
  v <- X %*% vcov(fit, regcoef.only=TRUE) %*% t(X)
  ndf <- if(is.matrix(v)) nrow(v) else 1
  se <- as.vector(if(ndf == 1) sqrt(v) else sqrt(diag(v)))
  Z <- est / se
  P <- if(length(idf)) 2 * pt(- abs(Z), idf) else 2 * pnorm(- abs(Z))
  if(conf.type != 'simultaneous') {
    if(length(bcoef)) {
      best <- t(matxv(X, bcoef, bmat=TRUE))
      lim <- bootBCa(est, best, type=boot.type, n=nobs(fit), seed=fit$seed,
                     conf.int=conf.int)
      if(is.matrix(lim)) {
        lower <- lim[1,]
        upper <- lim[2,]
      } else {
        lower <- lim[1]
        upper <- lim[2]
      }
    } else {
      lower <- est - zcrit*se
      upper <- est + zcrit*se
    }
  } else {
    u <- confint(multcomp::glht(fite, X,
                      df=if(length(idf)) idf else 0),
                 level=conf.int)$confint
    lower <- u[,'lwr']
    upper <- u[,'upr']
  }
  PP <- NULL; posterior.summary=''
  }
  res <- list(Contrast=est, SE=se,
              Lower=lower, Upper=upper,
              Z=Z, Pvalue=P, PP=PP,
              var=v, df.residual=idf,
              X=X, 
              cnames=if(type=='average') NULL else cnames,
              nvary=length(vary),
              conf.type=conf.type, conf.int=conf.int,
              posterior.summary=posterior.summary,
              cdraws = cdraws)
  if(type != 'average') res <- c(vary, res)
  
  r <- qr(v, tol=tol)
  nonred <- r$pivot[1:r$rank]   # non-redundant contrasts
  redundant <- (1:length(est)) %nin% nonred
  res$redundant <- redundant
  
  if(type=='joint') {
    est <- est[!redundant]
    v <- v[!redundant, !redundant, drop=FALSE]
    res$jointstat <- as.vector(est %*% solve(v, tol=tol) %*% est)
  }
  
  structure(res, class='contrast.rms')
}

print.contrast.rms <- function(x, X=FALSE, fun=function(u) u,
                               jointonly=FALSE, prob=0.95, ...)
{
  # See if a result of fun= on a Bayesian fit
  if('esta' %in% names(x)) {
    esta <- x$esta
    estb <- x$estb
    f <- function(x) {
      hpd <- HPDint(x, prob)
      r <- c(mean(x), median(x), hpd)
      names(r) <- c('Posterior Mean', 'Posterior Median',
                    paste(c('Lower', 'Upper'), prob, 'HPD'))
      r
    }
    cat('\nPosterior Summaries for First X Settings\n\n')
    print(t(apply(esta, 2, f)))
    cat('\nPosterior Summaries for Second X Settings\n\n')
    print(t(apply(estb, 2, f)))
    cat('\nPosterior Summaries of First - Second\n\n')
    print(t(apply(esta - estb, 2, f)))
    return(invisible())
    }
  
  edf <- x$df.residual
  sn <- if(length(edf)) 't' else 'Z'
  pn <- if(length(edf)) 'Pr(>|t|)' else 'Pr(>|z|)'
  w <- x[1 : (x$nvary + 7)]
  isn <- sapply(w, is.null)
  w <- w[! isn]
  
  if(length(w$Z))      w$Z      <- round(w$Z, 2)
  if(length(w$Pvalue)) w$Pvalue <- round(w$Pvalue, 4)
  if(length(w$PP))     w$PP     <- round(w$PP, 4)
  if(length(w$PP))     pn       <- 'Pr(Contrast>0)'
  
  no <- names(w)
  no[no=='SE'] <- 'S.E.'
  no[no=='Z']  <- sn
  no[no %in% c('Pvalue', 'PP')] <- pn
  
  cnames <- x$cnames
  if(! length(cnames))
    cnames <- if(x$nvary) rep('', length(x[[1]])) else
                          as.character(1 : length(x[[1]]))
  if(any(x$redundant)) cnames <- paste(ifelse(x$redundant, '*', ' '), cnames)
  w <- data.frame(w, row.names=paste(format(1:length(cnames)), cnames, sep=''))
  w$Contrast <- fun(w$Contrast)
  if(! all(1:10 == fun(1:10))) w$SE <- rep(NA, length(w$SE))
  w$Lower    <- fun(w$Lower)
  w$Upper    <- fun(w$Upper)

  # Assign modified names to w
  names(w) <- no

  # Print w
  if(!jointonly) {
    ## print(as.matrix(w), quote=FALSE)
    print(w, ...)
    if(any(x$redundant)) cat('\nRedundant contrasts are denoted by *\n')
  }
  
  jstat <- x$jointstat
  if(length(jstat)) {
    cat('\nJoint test for all contrasts=0:\n\n')
    ndf <- sum(!x$redundant)
    if(length(edf)) {
      Fstat <- jstat / ndf
      Pval <- 1 - pf(Fstat, ndf, edf)
      cat('F(', ndf, ',', edf, ')=', round(Fstat,3),', P=', round(Pval,4),
          '\n', sep='')
    } else {
      Pval <- 1 - pchisq(jstat, ndf)
      cat('Chi-square=', round(jstat, 2),' with ', ndf, ' d.f.  P=',
          round(Pval, 4),'\n', sep='')
    }
  }
  if(!jointonly && length(edf))cat('\nError d.f.=',edf,'\n')
  if(x$posterior.summary == '')
    cat('\nConfidence intervals are', x$conf.int, x$conf.type,
        'intervals\n')
  else {
    cat('\nIntervals are', x$conf.int, 'highest posterior density intervals\n')
    cat('Contrast is the posterior', x$posterior.summary, '\n')
    }
    
  if(X) {
    cat('\nDesign Matrix for Contrasts\n\n')
    if(is.matrix(x$X)) dimnames(x$X) <- list(cnames, dimnames(x$X)[[2]])
    print(x$X)
  }
  invisible()
}
