contrast <- function(fit, ...) UseMethod("contrast")

contrast.rms <-
  function(fit, a, b, a2, b2, ycut=NULL, cnames=NULL, fun=NULL, funint=TRUE,
           type=c('individual','average','joint'),
           conf.type=c('individual','simultaneous','profile'), usebootcoef=TRUE,
           boot.type=c('percentile','bca','basic'),
           posterior.summary=c('mean', 'median', 'mode'),
           weights='equal', conf.int=0.95, tol=1e-7, expand=TRUE, se_factor=4,
           plot_profile=FALSE, ...)
{
  type              <- match.arg(type)
  conf.type         <- match.arg(conf.type)
  boot.type         <- match.arg(boot.type)
  posterior.summary <- match.arg(posterior.summary)

  draws <- fit$draws
  bayes <- length(draws) > 0

  if(conf.type == 'profile' && type != 'individual')
    stop('conf.type=profile only works with type=individual')
  if(bayes & (type == 'joint' || conf.type == 'simultaneous'))
    stop('type=joint or conf.type=simultaneous not allowed for Bayesian models')
  
  zcrit <- if(length(idf <- fit$df.residual)) qt((1 + conf.int) / 2, idf) else
              qnorm((1 + conf.int) / 2)

  bcoef <- if(usebootcoef) fit$boot.Coef

  pmode <- function(x) {
    dens <- density(x)
    dens$x[which.max(dens$y)[1]]
    }
  
  if(! bayes) {
    betas <- coef(fit)
    iparm <- 1 : length(betas)
  }
  fite   <- fit
  ordfit <- inherits(fit, 'orm') || inherits(fit, 'lrm')
  if(ordfit) {
    nrp <- 1
    ## Note: is 1 for orm because vcov defaults to intercepts='mid' and
    ## we are overriding the default vcov uses for lrm
    fit$override_vcov_intercept <- 'mid'
    iparm <- c(fit$interceptRef, (num.intercepts(fit) + 1) : length(betas))
    betas <- betas[iparm]
    fite$coefficients <- betas    # for simult confint
    if(usebootcoef) bcoef <- bcoef[, iparm, drop=FALSE]
  } else nrp <- num.intercepts(fit, 'var')
    
  if(length(bcoef) && conf.type != 'simultaneous')
    conf.type <- switch(boot.type,
                        percentile = 'bootstrap nonparametric percentile',
                        bca        = 'bootstrap BCa',
                        basic      = 'basic bootstrap')

  partialpo <- inherits(fit, 'blrm') && fit$pppo > 0
  if(partialpo & ! length(ycut))
    stop('must specify ycut for partial prop. odds model')
  cppo      <- fit$cppo
  if(partialpo && ! length(cppo))
    stop('only implemented for constrained partial PO models')
  
  pred <- function(d) {
    ## predict.blrm duplicates rows of design matrix for partial PO models
    ## if ycut has length > 1 and only one observation is being predicted
    if(partialpo) predict(fit, d, type='x', ycut=ycut)
         else
           predict(fit, d, type='x')
    }
  
  da <- do.call('gendata', list(fit, factors=a, expand=expand))
  xa <- pred(da)
  if(! missing(b)) {
    db <- do.call('gendata', list(fit, factors=b, expand=expand))
    xb <- pred(db)
    }

  ma <- nrow(xa)

  if(missing(b)) {
    xb <- 0 * xa
    db <- da
  }
  mb <- nrow(xb)

  if(! missing(a2)) {
    if(missing(b) || missing(b2)) stop('b and b2 must be given if a2 is given')
    da2 <- do.call('gendata', list(fit, factors=a2, expand=expand))
    xa2 <- pred(da2)
    ma2 <- nrow(xa2)
    db2 <- do.call('gendata', list(fit, factors=b2, expand=expand))
    xb2 <- pred(db2)
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
          all(as.character(db[[w]]) == as.character(db2[[w]]))
      ## was da[[2]] 2023-09-07
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
    if(!missing(ycut)) stop('ycut not used with fun=')

    pa <- predict(fit, da, fun=fun, funint=funint, posterior.summary='all')
    pb <- predict(fit, db, fun=fun, funint=funint, posterior.summary='all')

    if(length(cnames)) colnames(pa) <- colnames(pb) <- cnames
    # If fun has an intercepts argument, the intecept vector must be
    # updated for each draw
    if(! length(cnames))
      cnames <- if(length(vary)) rep('', ncol(pa)) else
                          as.character(1 : ncol(pa))
    colnames(pa) <- colnames(pb) <- cnames

    res <- list(esta=pa, estb=pb,
                Xa=xa, Xb=xb, 
                nvary=length(vary))
    return(structure(res, class='contrast.rms'))
    }   # end if bayes & length(fun) ...
  
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
    ci     <- apply(cdraws, 2, rmsb::HPDint, prob=conf.int)
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
  ## vcov(lrm fit) has an override to use middle intercept - see above
  ndf <- if(is.matrix(v)) nrow(v) else 1
  se <- as.vector(if(ndf == 1) sqrt(v) else sqrt(Matrix::diag(v)))
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
    } else if(conf.type == 'profile') {
      w <- rms_profile_ci(X, fit, conf.int, est, se, plot_profile=plot_profile,
                          se_factor=se_factor, ...)
      lower <- w$lower
      upper <- w$upper
      LR    <- w$LR
      P     <- w$P
    } else {
      lower <- est - zcrit*se
      upper <- est + zcrit*se
    }
  } else {
    if(ordfit) {
      # glht uses vcov(fite) which for lrm & orm are sparse Matrix objects
      fite$non.slopes   <- 1L
      fite$interceptRef <- 1L
     if(! length(fite$var))
        fite$var <- Matrix::as.matrix(infoMxop(fite$info.matrix, i=iparm))
    }
    u <- confint(multcomp::glht(fite, X,
                      df=if(length(idf)) idf else 0),
                 level=conf.int)$confint
    lower <- u[,'lwr']
    upper <- u[,'upr']
  }
  PP <- NULL; posterior.summary=''
  }
  if(type != 'average' && length(ycut))
    cnames <- paste0(cnames, ' ', fit$yname, '=', ycut)
  res <- list(Contrast=est, SE=se,
              Lower=lower, Upper=upper,
              Z=Z, Pvalue=P, PP=PP,
              var=v, df.residual=idf,
              X=X, ycut=ycut, yname=fit$yname,  # was =if(length(ycut)) fit$yname
              cnames=if(type=='average') NULL else cnames,
              nvary=length(vary),
              conf.type=conf.type, conf.int=conf.int,
              posterior.summary=posterior.summary,
              cdraws = cdraws)
  if(conf.type == 'profile') res$LR <- LR
  if(type != 'average')      res <- c(vary, res)
  
  r <- qr(v, tol=tol)
  nonred <- r$pivot[1 : r$rank]   # non-redundant contrasts
  redundant <- (1 : length(est)) %nin% nonred
  res$redundant <- redundant
  
  if(type=='joint') {
    est <- est[! redundant]
    v <- v[! redundant, ! redundant, drop=FALSE]
    res$jointstat <- as.vector(est %*% solve(v, est, tol=tol))
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
      hpd <- rmsb::HPDint(x, prob)
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
  sn <- if(length(edf)) 't' else if(x$conf.type == 'profile') '\u03A7\u00B2' else 'Z'
  pn <- if(length(edf)) 'Pr(>|t|)' else if(x$conf.type == 'profile') 'Pr(>\u03A7\u00B2)' else 'Pr(>|z|)'
  if(length(x$LR)) {
    x$Z  <- x$LR
    x$LR <- NULL
  }
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
  if(length(x$y)) {
    w$.y. <- x$y
    names(w)[names(w) == '.y.'] <- x$yname
    }
  w$Contrast <- fun(w$Contrast)
  if(! all(1:10 == fun(1:10))) w$SE <- rep(NA, length(w$SE))
  w$Lower    <- fun(w$Lower)
  w$Upper    <- fun(w$Upper)

  # Assign modified names to w
  names(w) <- no

  if(x$conf.type == 'profile')             w$S.E. <- NULL
  if(length(w$S.E.) && all(is.na(w$S.E.))) w$S.E. <- NULL

  # Print w
  if(! jointonly) {
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
  if(!jointonly && length(edf)) cat('\nError d.f.=',edf,'\n')
  cotype <- if(x$conf.type == 'profile') 'profile likelihood' else x$conf.type
  if(x$posterior.summary == '')
    cat('\nConfidence intervals are', x$conf.int, cotype,
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

rms_profile_ci <-
  function(C, fit, conf.int, est_C, se_C, se_factor=4e0,
          plot_profile=FALSE, ...) {
  # Separate likelihood profile confidence intervals for contrasts in
  # each row of C.  est_C is estimated contrast, se_C is its standard error

  if(any(c('x', 'y') %nin% names(fit)))
    stop('to use profile likelihood you must specify x=TRUE, y=TRUE when fitting')
  X <- fit[['x']]

  crit  <- qchisq(conf.int, 1)
  p     <- ncol(C)
  m     <- nrow(C)

  if(p == (1 + length(fit$coefficients) - num.intercepts(fit))) C <- C[, -1, drop=FALSE]
  lower <- upper <- LR <- numeric(m)
  odev  <- getDeviance(fit)     # original deviance for full model
  odev  <- odev[length(odev)]

  g <- function(theta) {
    dev <- quickRefit(fit, X=Z[, -1, drop=FALSE], offset=theta * Z[, 1],
                      what='deviance', ...)
    if(is.list(dev) && length(dev$fail) && dev$fail) {
      message('Fit failed in profile likelihood.  theta=', format(theta), ' S.E.=', format(se),
           ' range of offsets:', paste(format(range(theta * Z[, 1])), collapse=', '))
      return(NA)
    }
    dev - odev - crit
  }

  p   <- ncol(C)

  for(i in 1 : m) {
    D <- C[i, , drop=FALSE]
    est <- est_C[i]
    se  <- se_C[i]
    v <- svd(rbind(D, diag(p)))$v
    u <- sqrt(sum(D ^ 2))
    beta_contrast <- v[, 1] * u
    # SVD has an arbitrary sign
    if(max(abs(D - beta_contrast)) > 1e-6) {
      v <- -v
      beta_contrast <- v[, 1] * u
      if(max(abs(D - beta_contrast)) > 1e-6)
        stop('SVD-generated contrast could not reproduce original contrast')
    }
    # Compute contrast to put on design matrix that gives the above contrast in betas
    v <- v / u
    Z <- X %*% v
    # Likelihood ratio chi-square obtained by removing first column of Z
    drop1 <- quickRefit(fit, X=Z[, -1, drop=FALSE], what='deviance', ...)
    drop1 <- drop1[length(drop1)]
    LR[i] <- drop1 - odev
    if(plot_profile) {
      thetas      <- seq(est - se_factor * se, est + se_factor * se, length=50)
      ch_deviance <- rep(NA, length(thetas))
      for(j in 1 : length(thetas)) ch_deviance[j] <- g(thetas[j])
      plot(thetas, ch_deviance, xlab='Contrast Estimate',
           ylab='Change in Deviance From Full Model')
      abline(v=c(est - se, est, est + se), col='blue')
      title(paste('Contrast', i))
      title(sub='Vertical lines are at point estimate of contrast \u00b1 S.E.', adj=1, cex.sub=0.65)  
    }
    hi <- try(uniroot(g, c(est + se/100, est + se_factor * se))$root)
    if(inherits(hi, 'try-error')) hi <-  Inf
    lo <- try(uniroot(g, c(est - se_factor * se, est - se/100))$root)
    if(inherits(lo, 'try-error')) lo <- -Inf
    lower[i] <- lo
    upper[i] <- hi
  }
list(lower=lower, upper=upper, LR=LR, P=1. - pchisq(LR, 1))
}

