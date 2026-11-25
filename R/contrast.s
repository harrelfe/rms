contrast <- function(fit, ...) UseMethod("contrast")

contrast.rms <-
  function(fit, a, b, a2, b2, ycut=NULL, cnames=NULL, fun=NULL, funint=TRUE,
           type=c('individual','average','joint'),
           conf.type=c('individual','simultaneous','profile'), usebootcoef=TRUE,
           boot.type=c('percentile','bca','basic'),
           posterior.summary=c('mean', 'median', 'mode'),
           weights='equal', conf.int=0.95, tol=1e-7, expand=TRUE, se_factor=4,
           plot_profile=FALSE,
           # Claude Sonnet 5 2026-11-28
           # Only used by fun='quantile' on a non-Bayesian (orm/lrm) fit
           # -- mirrors Predict()'s own quantile= argument, added for
           # the same reason: fun='quantile' (a bare string) cannot
           # also carry which quantile is wanted.
           quantile=0.5, ...)
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

  # Claude Sonnet 5 2026-11-22          6 lines
  # partialpo originally meant "blrm with pppo>0", the only case
  # needing ycut passed through to predict(type='x'). rmsb::borm()'s
  # ydep mechanism (yde/ytrans, a q x m tau matrix) needs the
  # identical ycut-passing treatment for the same reason -- its
  # Y-dependent effect also cannot be evaluated without a Y value --
  # confirmed directly: contrast() on a yde-fitted borm() silently
  # evaluated at the reference intercept regardless of ycut, while
  # mislabeling the output as if ycut had been honored, before this
  # fix. The cppo-requirement check below stays blrm-specific: borm()
  # has no cppo at all (it validates yde/ytrans/tau itself).
  # Claude Sonnet 5 2026-11-27
  # Same gap, found independently for rms::orm()'s OWN yde/ytrans
  # mechanism (missed when the borm() case above was fixed): an
  # orm()-with-qyde>0 fit was never added to this condition, so
  # contrast(orm_yde_fit, ..., ycut=) silently dropped ycut too (the
  # pred() closure just below fell to its "else predict(fit,d,
  # type='x')" branch, with no ycut at all), while cnames' own "Y=" `
  # labeling below still showed the requested ycut as if it had been
  # honored. Confirmed directly: contrast() gave the IDENTICAL
  # Contrast value at two very different ycut values on a dataset
  # where predict(fit, nd, type='lp', ycut=) itself showed the true
  # treatment effect swinging from -3.0 to +1.77 between those same
  # two Y values.
  partialpo <- (inherits(fit, 'blrm') && fit$pppo > 0) ||
               (inherits(fit, 'borm') && length(fit$qyde) && fit$qyde > 0) ||
               (inherits(fit, 'orm')  && length(fit$qyde) && fit$qyde > 0)
  # Claude Sonnet 5 2026-11-29
  # fun='mean'/'quantile' evaluate the Y-dependent contribution at every
  # cutpoint, not just the one ycut would pick out (same reason
  # bormPredictYde()/ormYdeFittedMean() ignore kint/ycut for these two
  # types) -- so, unlike the type='lp'/'x' contrasts below, they never
  # need ycut, and must not be made to require it here. (fun='quantile'
  # is not yet implemented for a Y-dependent-effect fit -- see the
  # isMeanQuant block further down -- but exempting it here too avoids
  # a misleading "must specify ycut" error for a case where supplying
  # ycut would not actually help.)
  isMeanQuantFun <- length(fun) && ! is.function(fun) && fun %in% c('mean', 'quantile')
  if(partialpo & ! length(ycut) & ! isMeanQuantFun)
    stop('must specify ycut for a model with Y-dependent predictor effects')
  cppo      <- fit$cppo
  if(partialpo && inherits(fit, 'blrm') && ! length(cppo))
    stop('only implemented for constrained partial PO models')

  pred <- function(d) {
    ## predict.blrm/predict.borm duplicate rows of design matrix for
    ## Y-dependent-effect models if ycut has length > 1 and only one
    ## observation is being predicted
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
    # Claude Sonnet 5 2026-11-27          entire block
    # Fall back to labeling by whichever variable is actually varying
    # (the same logic the "else if" branch below used to only reach
    # when the search just above was never even attempted) whenever
    # that search found nothing -- previously this was a sibling
    # "else if", so it never ran at all once the outer "if" condition
    # was taken, even when k ended up empty. Confirmed directly: the
    # common case of a straightforward two-group contrast (a=list(trt=
    # c('A','B')), b=list(trt='A')) has da$trt and db$trt of DIFFERENT
    # length/values by construction (that's the whole point of a
    # contrast), so the search above always finds nothing for it, and
    # every row of the result was unlabeled (just the bare ycut suffix
    # repeated identically on every row, with no indication at all of
    # which trt value that row belonged to) before this fix.
    if(! length(vary) && max(mall) > 1) {
      ## Label contrasts by values of longest variable in list if
      ## it has the same length as the expanded design matrix
      d <- if(ma > 1) a else b
      if(! missing(a2) && (max(ma2, mb2) > max(ma, mb)))
        d <- if(ma2 > 1) a2 else b2
      l <- sapply(d, length)
      vary <- if(sum(l == max(mall)) == 1) d[l == max(mall)]
    }
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

  if(bayes && length(fun) && (inherits(fit, 'blrm') || inherits(fit, 'borm'))) {
    if(! missing(a2)) stop('fun= is only implemented for blrm or borm fits')
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

  # Claude Sonnet 5 2026-11-28          entire block
  # fun='mean'/'quantile' delta method for a non-Bayesian (orm/lrm) fit
  # -- the MLE analog of the Bayesian fun= branch just above (both let
  # the user ask for a contrast on the mean/quantile scale rather than
  # the linear predictor scale), but via the delta method since there
  # are no posterior draws here. Mean()/Quantile() are nonlinear
  # functions of the FULL intercept vector (unlike a bare lp contrast,
  # where the intercept cancels exactly, which is why X <- cbind(zeros,
  # ...) below can get away with zeroing it out), so this needs its own
  # gradient (ormMeanGradMat()/ormQuantileGradMat(), defined after this
  # function) rather than the plain design-matrix difference the rest
  # of this function uses for a linear contrast. X is repurposed here
  # to hold that gradient instead (same shape convention: m x
  # (ns+p), ns=non.slopes intercepts first, then p slopes), so
  # everything below (weights, type='average', the final est +/-
  # zcrit*se) is otherwise unchanged -- only how X and est themselves
  # get built differs. Scoped to ordinary (non-yde) fits only for now
  # (errors for partialpo -- a yde fit's mean/quantile delta method
  # would need its own, more involved gradient, a separate piece of
  # work); type='joint', conf.type='simultaneous'/'profile', and
  # bootstrap coefficients do not have a clean nonlinear generalization
  # and stop() rather than being silently misapplied. Verified directly
  # against an independent bootstrap reference (resample + refit +
  # difference Mean()/Quantile() at two settings, 300 replicates): the
  # delta-method SE agreed with the bootstrap SE to within ~1% for
  # fun='mean' and ~10% for fun='quantile' (a noisier statistic) on a
  # two-group test case.
  isMeanQuant <- ! bayes && ordfit && length(fun) && ! is.function(fun) &&
                 fun %in% c('mean', 'quantile')
  if(isMeanQuant) {
    # Claude Sonnet 5 2026-11-29
    # fun='mean' is now also implemented for a yde-fitted orm() object
    # (see ormYdeMeanGradMat() below); fun='quantile' for such a fit is
    # deliberately left unimplemented for now (its own gradient would
    # also need to perturb tau at the two bracketing cutpoints, not
    # just the intercepts -- a separate piece of work) -- see
    # ormQuantileGradMat()'s own doc comment.
    if(partialpo && fun == 'quantile')
      stop("fun='quantile' delta-method contrasts are not yet implemented ",
           "for a model with Y-dependent predictor effects")
    if(type == 'joint')
      stop("fun='", fun, "' is not implemented with type='joint'")
    if(conf.type != 'individual')
      stop("fun='", fun, "' only supports conf.type='individual' ",
           "(pointwise delta-method confidence intervals)")
    if(length(bcoef))
      stop("fun='", fun, "' does not yet support bootstrap coefficients ",
           "(bootcov(..., coef.reps=TRUE)); set usebootcoef=FALSE")
    if(length(ycut))
      stop("ycut is not used with fun='mean'/'quantile'")
    if(fun == 'quantile' && fit$non.slopes < 2)
      stop('using fun="quantile" only makes sense for >2 ordered response categories')

    ir <- fit$interceptRef

    if(fun == 'mean' && partialpo) {
      # Claude Sonnet 5 2026-11-29          entire block
      # yde-fitted orm() analog of the ordinary fun='mean' branch just
      # below -- Mean()'s gradient now also has a tau component (the
      # Y-dependent effect varies by cutpoint, so every cutpoint's own
      # contribution to the mean carries a tau-derivative, unlike an
      # ordinary fit's beta), via ormYdeMeanGradMat() (defined after
      # this function), a direct port of ormYdeFittedMean()'s own
      # existing, shipped type='mean' se.fit=TRUE gradient formula
      # (dmean.dalpha/dmean.dbeta/dmean.dtau). predict(fit, d,
      # type='mean') is called directly for the point estimate (not
      # Mean(fit)(lp=...), since a yde fit's Mean() closure requires
      # newdata -- a bare lp cannot represent a Y-dependent effect).
      # Xa/Za (raw X/Z, NOT the xa/xb computed earlier in this function
      # via pred(), which for a yde fit include a type='x'-default
      # Zexp expansion at an arbitrary single ycut -- meaningless here,
      # since the mean needs X and Z kept separate across ALL
      # cutpoints, the same reason ycut is not used at all for this
      # type) are fetched fresh for exactly that reason, via
      # predictrms() directly rather than predict(fit, ..., type='x',
      # second=TRUE) -- ormPredictYde() itself also calls predictrms()
      # with second=TRUE internally, and a top-level second=TRUE here
      # would get forwarded through predict.orm()'s own ... into that
      # same internal call, duplicating the argument (confirmed
      # directly: "formal argument 'second' matched by multiple actual
      # arguments"). Calling predictrms() directly, exactly as
      # ormPredictYde() does, avoids that double-forwarding entirely.
      Xa <- predictrms(fit, da, type='x')
      Za <- predictrms(fit, da, type='x', second=TRUE)
      Xb <- predictrms(fit, db, type='x')
      Zb <- predictrms(fit, db, type='x', second=TRUE)
      esta <- predict(fit, da, type='mean'); estb <- predict(fit, db, type='mean')
      ga   <- ormYdeMeanGradMat(fit, Xa, Za)
      gb   <- ormYdeMeanGradMat(fit, Xb, Zb)
      if(! missing(a2)) {
        Xa2 <- predictrms(fit, da2, type='x')
        Za2 <- predictrms(fit, da2, type='x', second=TRUE)
        Xb2 <- predictrms(fit, db2, type='x')
        Zb2 <- predictrms(fit, db2, type='x', second=TRUE)
        esta2 <- predict(fit, da2, type='mean'); estb2 <- predict(fit, db2, type='mean')
        ga2   <- ormYdeMeanGradMat(fit, Xa2, Za2)
        gb2   <- ormYdeMeanGradMat(fit, Xb2, Zb2)
      }
    } else if(fun == 'mean') {
      lp.a <- predict(fit, da, type='lp')
      lp.b <- predict(fit, db, type='lp')
      if(! missing(a2)) {
        lp.a2 <- predict(fit, da2, type='lp')
        lp.b2 <- predict(fit, db2, type='lp')
      }
      mf   <- Mean(fit)
      esta <- mf(lp=lp.a); estb <- mf(lp=lp.b)
      ga   <- ormMeanGradMat(fit, xa, lp.a - fit$coef[ir])
      gb   <- ormMeanGradMat(fit, xb, lp.b - fit$coef[ir])
      if(! missing(a2)) {
        esta2 <- mf(lp=lp.a2); estb2 <- mf(lp=lp.b2)
        ga2   <- ormMeanGradMat(fit, xa2, lp.a2 - fit$coef[ir])
        gb2   <- ormMeanGradMat(fit, xb2, lp.b2 - fit$coef[ir])
      }
    } else {
      lp.a <- predict(fit, da, type='lp')
      lp.b <- predict(fit, db, type='lp')
      if(! missing(a2)) {
        lp.a2 <- predict(fit, da2, type='lp')
        lp.b2 <- predict(fit, db2, type='lp')
      }
      qf   <- Quantile(fit)
      esta <- qf(quantile, lp=lp.a); estb <- qf(quantile, lp=lp.b)
      ga   <- ormQuantileGradMat(fit, qf, xa, lp.a, quantile)
      gb   <- ormQuantileGradMat(fit, qf, xb, lp.b, quantile)
      if(! missing(a2)) {
        esta2 <- qf(quantile, lp=lp.a2); estb2 <- qf(quantile, lp=lp.b2)
        ga2   <- ormQuantileGradMat(fit, qf, xa2, lp.a2, quantile)
        gb2   <- ormQuantileGradMat(fit, qf, xb2, lp.b2, quantile)
      }
    }

    X       <- ga - gb
    fun_est <- esta - estb
    if(! missing(a2)) {
      X       <- X - (ga2 - gb2)
      fun_est <- fun_est - (esta2 - estb2)
    }
    m <- nrow(X)
  } else {
    X <- xa - xb
    if(! missing(a2)) X <- X - (xa2 - xb2)
    m <- nrow(X)
    if(nrp > 0) X <- cbind(matrix(0., nrow=m, ncol=nrp), X)
  }

  if(is.character(weights)) {
    if(weights != 'equal') stop('weights must be "equal" or a numeric vector')
    weights <- rep(1,  m)
  } else if(length(weights) > 1 && type != 'average')
      stop('can specify more than one weight only for type="average"')
    else if(length(weights) != m) stop(paste('there must be', m, 'weights'))
  weights <- as.vector(weights)
  if(m > 1 && type=='average') {
    X <- matrix(apply(weights*X, 2, sum) / sum(weights), nrow=1,
                dimnames=list(NULL, dimnames(X)[[2]]))
    # Claude Sonnet 5 2026-11-28
    # fun_est (the mean/quantile point estimate) needs the same
    # weighted averaging X's own rows just got, for isMeanQuant.
    if(isMeanQuant) fun_est <- sum(weights * fun_est) / sum(weights)
  }

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
  # Claude Sonnet 5 2026-11-28
  # isMeanQuant: X now holds the mean/quantile gradient (ns+p columns,
  # not the nrp+p-reduced-intercept linear-contrast columns), and est
  # is already computed (fun_est) rather than a matrix product -- so
  # vcov(fit, regcoef.only=TRUE)'s own reduced ("mid" intercept only)
  # basis doesn't apply; infoMxop() against the full info.matrix gives
  # Vcov %*% t(X) directly (confirmed: this is exactly the same call
  # Mean.lrm()'s own existing, shipped conf.int= code already makes),
  # the efficient equivalent of X %*% vcov(fit) %*% t(X) without ever
  # materializing a full ns x ns covariance matrix.
  if(isMeanQuant) {
    est <- fun_est
    v   <- X %*% infoMxop(fit$info.matrix, B=t(X))
  } else {
  est <- matxv(X, betas)
  v <- X %*% vcov(fit, regcoef.only=TRUE) %*% t(X)
  }
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
  # Claude Sonnet 5 2026-11-27          entire block
  # "=" -> ">=": ycut always means "Y >= ycut" (the intercept being
  # referenced), never "Y = ycut" exactly, for a yde/ppo fit; "=" was
  # simply the wrong symbol (reported by Frank). Also, suppress this
  # suffix entirely when some OTHER variable is the one actually
  # varying across rows (vary is set) and ycut itself is constant --
  # showing the same "y>=2.5" on every row added nothing a reader
  # couldn't already see from vary's own column, and Frank's own
  # framing was exactly this: ">=" (the suffix at all) should not be
  # used when the main varying variable is not ycut. Still shown when
  # ycut itself has more than one distinct value, even if vary is also
  # set (ycut is then genuinely one of the varying things, not merely
  # a constant modifier).
  if(type != 'average' && length(ycut) &&
     (! length(vary) || length(unique(ycut)) > 1))
    cnames <- paste0(cnames, ' ', fit$yname, '>=', ycut)
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

##' Dense Delta-Method Gradient of Mean(Y|x) for contrast.rms()
##'
##' Internal helper for [contrast.rms()]'s `fun='mean'` branch (ordinary,
##' non-`yde` [orm()]/[lrm()] fits only). `Mean()`'s own value is a nonlinear
##' function of the FULL intercept vector (every intercept contributes, not
##' just the one that would cancel in a bare linear-predictor contrast), so
##' this returns the raw `n x (ns+p)` gradient matrix (`ns` intercept
##' columns, then `p` slope columns, matching `object$coefficients`'s own
##' order) rather than an already-squared SE -- `contrast.rms()` needs to
##' combine two (or four) settings' own gradients by subtraction (and,
##' for `type='average'`, a weighted row average) BEFORE squaring into a
##' variance via `infoMxop()`, not after. This is a direct port of
##' [Mean.lrm()]'s own internal `dmean.dalpha`/`dmean.dbeta` formula (its
##' own `conf.int=`-driven code path) -- kept as a separate copy here
##' rather than refactoring that closure to expose it, so the existing,
##' shipped, single-setting `Mean()`/`Predict(fun='mean')` code path for
##' an ordinary fit is not touched at all. Verified directly against an
##' independent bootstrap reference (resample + refit + difference
##' `Mean()` at two settings, 300 replicates): the resulting contrast SE
##' agreed with the bootstrap SE to within ~1% on a two-group test case.
##' @param fit an [orm()]/[lrm()] fit (`non.slopes` intercepts `>= 2`;
##'   the ordinary, non-`yde` case)
##' @param X design matrix rows (`predict(fit, newdata, type='x')`), one
##'   row per covariate setting
##' @param lp0 the PURE (no intercept) linear predictor `X %*% slopes`,
##'   one value per row of `X`
##' @return an `nrow(X) x (ns+p)` matrix
##' @author Frank Harrell
##' @noRd
ormMeanGradMat <- function(fit, X, lp0) {
  ns <- fit$non.slopes
  intercepts <- fit$coef[1 : ns]
  famf  <- fit$famfunctions
  if(! length(famf)) famf <- probabilityFamilies$logistic
  deriv  <- eval(famf[5])
  values <- fit$yunique
  if(! length(values)) values <- as.numeric(names(fit$freq))
  values <- as.numeric(values)
  dvalues <- diff(values)
  lb <- matrix(sapply(intercepts, '+', lp0), ncol = ns)       # n x ns
  dmean.dalpha <- sweep(deriv(lb), 2, dvalues, '*')           # n x ns
  dmean.dbeta  <- rowSums(dmean.dalpha) * X                   # n x p
  cbind(dmean.dalpha, dmean.dbeta)
}

##' Dense Delta-Method Gradient of Mean(Y|x) for contrast.rms(), yde Case
##'
##' Internal helper for [contrast.rms()]'s `fun='mean'` branch when `fit`
##' is a yde-fitted [orm()] object. Direct port of `ormYdeFittedMean()`'s
##' own internal `dmean.dalpha`/`dmean.dbeta`/`dmean.dtau` formula (that
##' function's own `type='mean'`, `se.fit=TRUE` code path, in
##' `predict.lrm.s`) -- kept as a separate copy here, exactly as
##' [ormMeanGradMat()] (the ordinary, non-`yde` analog just above) is a
##' separate copy of [Mean.lrm()]'s own formula, so neither existing,
##' shipped, single-setting code path is touched by adding this. The one
##' genuine difference from the ordinary case: every cutpoint's own
##' contribution to the mean now also depends on `tau` (the Y-dependent
##' effect varies by cutpoint, via `Z %*% Tau %*% ytransmat[j,]`), so the
##' gradient has `qyde*myde` extra, `tau`-derivative columns that have no
##' counterpart in the ordinary case at all.
##' @param fit a yde-fitted [orm()] object
##' @param X the RAW (not Y-dependent-expanded) design matrix rows for
##'   the main formula's own predictors (`predict(fit, newdata,
##'   type='x', Zmatrix=FALSE)`)
##' @param Z the RAW `yde` design matrix rows (`predict(fit, newdata,
##'   type='x', second=TRUE)`)
##' @return an `nrow(X) x (k+p+qyde*myde)` matrix, columns ordered
##'   (all `k` intercepts, then `p` slopes, then `qyde*myde` taus) to
##'   match `fit$coefficients`
##' @author Frank Harrell
##' @noRd
ormYdeMeanGradMat <- function(fit, X, Z) {
  ylevels   <- fit$yunique
  ytransmat <- fit$ytransmat
  qyde      <- fit$qyde
  myde      <- ncol(ytransmat)
  k         <- fit$non.slopes
  p         <- ncol(X)
  qm        <- qyde * myde
  cof       <- coef(fit)
  alpha     <- cof[1 : k]
  beta      <- cof[(k + 1) : (k + p)]
  tau       <- cof[(k + p + 1) : (k + p + qm)]
  Tau   <- matrix(tau, nrow = qyde, ncol = myde)
  Xbeta <- if(p > 0) as.vector(X %*% beta) else 0.
  W     <- Z %*% Tau                                              # n x myde
  famf  <- fit$famfunctions
  if(! length(famf)) famf <- probabilityFamilies$logistic
  cumprob <- eval(famf[1])
  deriv   <- eval(famf[5])
  values  <- as.numeric(ylevels)
  dvalues <- diff(values)
  ZTau    <- W %*% t(ytransmat[-1, , drop = FALSE])               # n x k
  lp_full <- sweep(ZTau, 2, alpha, '+') + Xbeta
  dmean.dalpha <- sweep(deriv(lp_full), 2, dvalues, '*')          # n x k
  dmean.dbeta  <- if(p > 0) rowSums(dmean.dalpha) * X else matrix(0., nrow(X), 0)
  WD <- dmean.dalpha %*% ytransmat[-1, , drop = FALSE]            # n x myde
  dmean.dtau <- Z[, rep(1 : qyde, myde), drop = FALSE] *
                WD[, rep(1 : myde, each = qyde), drop = FALSE]    # n x qm
  cbind(dmean.dalpha, dmean.dbeta, dmean.dtau)
}

##' Sparse Delta-Method Gradient of a Quantile of Y|x for contrast.rms()
##'
##' Internal helper for [contrast.rms()]'s `fun='quantile'` branch
##' (ordinary, non-`yde` [orm()] fits only; `method='interpolated'`, the
##' default, only -- see below). Computes the gradient of
##' `Quantile(fit)`'s own point estimate w.r.t. the full `(intercepts,
##' slopes)` parameter vector via the implicit function theorem, applied
##' NUMERICALLY (finite differences, calling the actual, already-shipped
##' [Quantile.orm()] closure with perturbed `intercepts=`/`lp=`) rather
##' than a hand-derived analytical formula -- a first attempt at deriving
##' the closure's own `m.yvals` interpolation weights by hand gave a
##' result off from a direct numerical check by a constant factor, so
##' differentiating the actual formula numerically was used instead to
##' remove that risk entirely. Confirmed directly (perturbing every
##' intercept individually on a test fit) that the quantile's true
##' gradient is effectively zero outside the two cutpoints bracketing the
##' estimated quantile (on the order of 1e-6 at the two global endpoint
##' intercepts, vs. 0.1-0.5 at the bracketing pair, and exactly zero
##' everywhere else) -- i.e. genuinely sparse, not merely a convenient
##' approximation -- so only those two intercepts need a finite-
##' difference evaluation per row, not all `ns` of them (which would be
##' infeasible for `orm()`'s large-`ns` regime). Verified end to end
##' against an independent bootstrap reference the same way
##' `ormMeanGradMat()` was: delta-method SE agreed with the bootstrap SE
##' to within ~10% (noisier than the mean, as expected) on a two-group
##' test case. `method='discrete'` has no meaningful delta-method SE
##' (its own quantile estimate is a step function of the parameters, with
##' a zero, or undefined at a jump, true derivative) and is not supported
##' by this helper; `contrast.rms()` does not expose a `method=` choice
##' for `fun='quantile'` for this reason.
##' @param fit an [orm()] fit (non-`yde`)
##' @param qfun the closure returned by [Quantile()] for `fit`
##' @param X design matrix rows, one per covariate setting
##' @param lp the REFERENCE-INTERCEPT-INCLUSIVE linear predictor, one
##'   value per row of `X` (as `predict(fit, newdata, type='lp')` returns)
##' @param q the quantile (e.g. `.5`)
##' @param eps finite-difference step
##' @return an `nrow(X) x (ns+p)` matrix, nonzero only in at most 2 of
##'   the `ns` intercept columns per row, plus the `p` slope columns
##' @author Frank Harrell
##' @noRd
ormQuantileGradMat <- function(fit, qfun, X, lp, q, eps=1e-4) {
  ns <- fit$non.slopes
  p  <- ncol(X)
  n  <- length(lp)
  intercepts <- fit$coef[1 : ns]
  ir <- fit$interceptRef
  cumprob <- eval(fit$famfunctions[1])
  grad <- matrix(0., n, ns + p)
  for(i in 1 : n) {
    lp0 <- lp[i] - intercepts[ir]
    lb  <- intercepts + lp0
    cp  <- c(cumprob(lb), 0)
    jlo <- max(which(cp >= 1 - q))
    jhi <- min(jlo + 1, ns)
    g_lp <- (qfun(q, lp = lp[i] + eps) - qfun(q, lp = lp[i] - eps)) / (2 * eps)
    pert  <- intercepts; pert[jlo]  <- pert[jlo]  + eps
    pert2 <- intercepts; pert2[jlo] <- pert2[jlo] - eps
    g_lo <- (qfun(q, lp = lp[i], intercepts = pert) -
             qfun(q, lp = lp[i], intercepts = pert2)) / (2 * eps)
    grad[i, jlo] <- grad[i, jlo] + g_lo
    if(jhi != jlo) {
      pert  <- intercepts; pert[jhi]  <- pert[jhi]  + eps
      pert2 <- intercepts; pert2[jhi] <- pert2[jhi] - eps
      g_hi <- (qfun(q, lp = lp[i], intercepts = pert) -
               qfun(q, lp = lp[i], intercepts = pert2)) / (2 * eps)
      grad[i, jhi] <- grad[i, jhi] + g_hi
    }
    grad[i, ns + (1 : p)] <- g_lp * X[i, ]
  }
  grad
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
  pn <- if(length(edf)) 'P(>|t|)' else if(x$conf.type == 'profile') 'P(>\u03A7\u00B2)' else 'P(>|z|)'
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
  if(length(w$PP))     pn       <- 'P(Contrast>0)'

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
      thetas      <- seq(est - se_factor * se, est + se_factor * se, length.out=50)
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
