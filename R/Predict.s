Predict <-
  function(object, ...,
           fun=NULL, funint=TRUE, type=c("predictions","model.frame","x"),
           np=200,
           conf.int=.95,
           conf.type=c('mean', 'individual', 'simultaneous'),
           usebootcoef=TRUE, boot.type=c('percentile', 'bca', 'basic'),
           posterior.summary=c('mean', 'median', 'mode'),
           adj.zero=FALSE, ref.zero=FALSE,
           kint=NULL, ycut=NULL, time=NULL,
           loglog=FALSE, digits=4, name, factors=NULL,
           offset=NULL, quantile=0.5)
{
  fit       <- object
  type      <- match.arg(type)
  conf.type <- match.arg(conf.type)
  boot.type <- match.arg(boot.type)
  posterior.summary <- match.arg(posterior.summary)
  draws     <- fit$draws
  bayes     <- length(draws) > 0
  if(bayes && conf.type == 'simultaneous')
    stop('conf.type simultaneous does not work for Bayesian models')
  
  oldopt <- options('digits')
  options(digits=digits)
  on.exit(options(oldopt))

  cl <- class(fit)
  isblrm <- 'blrm' %in% cl
  isborm <- 'borm' %in% cl
  isorm  <- 'orm'  %in% cl
  islrm  <- 'lrm'  %in% cl

  Center <- 0.
  if('cph' %in% cl) Center <- fit$center

  kintgiven <- length(kint) > 0
  if(! length(kint)) {
    kint <- fit$interceptRef
    if(! length(kint)) kint <- 1
  }

  # Claude Sonnet 5 2026-11-17          entire block
  # Generalized to also recognize a yde-fitted borm() object, which
  # has no pppo/cppo fields at all (borm()'s own non-PO mechanism is
  # qyde/ytrans/tau) -- previously this unconditionally forced pppo to
  # 0 for any non-blrm object, so a yde-fitted borm() silently lost
  # its Y-dependent effect here. blrm's own cppo-presence requirement
  # (only the constrained, cppo-based partial-PO case is supported by
  # this function) stays blrm-specific, since borm() has no
  # unconstrained-PO alternative to rule out -- qyde > 0 always means
  # the ytrans-constrained case for borm().
  pppo <- fit$pppo
  if(! isblrm) pppo <- 0
  qyde <- fit$qyde
  # Claude Sonnet 5 2026-11-17          1 line
  # Also recognize a yde-fitted orm() object (qyde > 0), not just borm().
  if((! isborm && ! isorm) || ! length(qyde)) qyde <- 0
  partialpo <- pppo > 0L || qyde > 0L
  if(isblrm) cppo <- fit$cppo
  if(partialpo && isblrm && ! length(cppo))
    stop('only implemented for constrained partial PO models')
  ylevels   <- if(isblrm || isborm) fit$ylevels else fit$yunique
  if(islrm || isorm || isblrm || isborm) {
    if(kintgiven && ! length(ycut)) ycut <- ylevels[kint + 1]
    if(length(ycut) && ! kintgiven)
      kint <- if(all.is.numeric(ylevels)) which(ylevels == ycut) - 1
              else max((1 : length(ylevels))[ylevels <= ycut]) - 1
  }

  # Claude Sonnet 5 2026-11-17          1 line
  # A yde-fitted orm() object routes through predict() (-> predict.orm()
  # -> ormPredictYde()), same as blrm()/borm(), not predictrms()
  # directly -- predictrms() alone knows nothing about tau/Z.
  isormyde  <- isorm  && qyde > 0L
  # Claude Sonnet 5 2026-11-22
  # Same reason isormyde needs its own branch below: bormPredictYde()'s
  # own type='lp' stops on fun= (a per-draw-then-summarize concept it
  # doesn't implement), so a yde-fitted borm() object needs the same
  # no-fun-passthrough treatment isormyde already gets, not the
  # isblrm/isborm branch's fun=fun forwarding.
  isbormyde <- isborm && qyde > 0L
  Pred <- function(type='lp') {
    if(type == 'x') {
      if(isblrm || isborm || isormyde) predict(fit, settings, type='x', kint=kint, ycut=ycut,
                         zcppo=FALSE)
      else
        predictrms(fit, settings, type='x')
    }
    else {
      # Claude Sonnet 5 2026-11-20          6 lines
      # isormyde/isbormyde split out of the blrm/borm branch: fun=/
      # funint=/posterior.summary= are a Bayesian-draws concept (apply
      # fun per draw, then summarize) that predict.orm()/ormPredictYde()
      # (an MLE fit, no draws) does not implement at all, and that
      # bormPredictYde() doesn't implement for fun= specifically (though
      # it does support posterior.summary=/cint= already) -- fun is
      # instead applied later, below, exactly as for an ordinary
      # (non-yde) orm/lrm fit.
      # Claude Sonnet 5 2026-11-22          entire block
      # fun=Mean(fit) (the real function), not fun=fun (the literal
      # string 'mean'), when isMean: predict.blrm()'s own 'lp' code has
      # no is.function() check before trying to call fun() in its
      # per-draw loop (confirmed directly, on a PLAIN non-yde blrm()
      # fit, unrelated to any yde/borm work: "could not find function
      # 'fun'"). Mean.blrm()'s closure signature
      # (lp, lptau=, intercepts=, codes=) is exactly what that loop
      # calls it with, so this isn't a workaround -- it's the real
      # function the string 'mean' is shorthand for, confirmed by
      # reading how that branch finishes: it HPD-summarizes the
      # already-fun-transformed per-draw values, so lower/upper come
      # back on the mean (not log-odds) scale automatically. Passing
      # NULL instead (an earlier, incomplete fix attempt) avoided the
      # crash but left lower/upper on the untransformed log-odds scale
      # while xb got transformed later via m(xb) -- confirmed directly
      # to produce a nonsensical result (the point estimate landing
      # outside its own "confidence interval").
      if(isblrm || (isborm && ! isbormyde)) predict(fit, settings, type='lp',
                       fun=if(isMean) Mean(fit) else fun, funint=funint,
                       kint=kint, ycut=ycut,
                       posterior.summary=posterior.summary, cint=conf.int)
      else if(isormyde) predict(fit, settings, type='lp', kint=kint, ycut=ycut,
                                 cint=conf.int)
      else if(isbormyde) predict(fit, settings, type='lp', kint=kint, ycut=ycut,
                                  posterior.summary=posterior.summary, cint=conf.int)
      else predictrms(fit, settings, kint=kint,
                      ref.zero=ref.zero,
                      type='lp', conf.int=conf.int, conf.type=conf.type)
    }
  }
  
  dotlist <- if(length(factors)) factors else rmsArgs(substitute(list(...)))
  fname   <- if(missing(name)) '' else name
  at      <- fit$Design
  assume  <- at$assume.code
  name    <- at$name	##interactions are placed at end by design

  ioff <- attr(fit$terms, 'offset')
  if(length(ioff)) {
    offsetExpression <- rownames(attr(fit$terms, 'factors'))[ioff]
    offsetVariableName <- all.vars(as.formula(paste('~', offsetExpression)))
    if(! length(offset))
      stop('model has offset term but offset=list(...) not given to Predict')
    if(length(offset) > 1) stop('offset may only contain one variable')
    if(length(offset[[1]]) != 1) stop('offset variable must contain 1 value')
    if(names(offset) != offsetVariableName)
      stop(paste('offset does not have correct variable name (',
                 offsetVariableName, ')', sep=''))
  }
    
  if('time' %in% name) {
    dotlist$time <- time
    time <- NULL
  }

  if(length(fname) > 1 || (length(dotlist) == 0 && fname == '')) {
    m <- match.call(expand.dots=FALSE)
    m[[1]] <- as.name('Predict')
    nams <- if(length(fname) > 1) fname else name[assume != 9]
    res <- vector('list', length(nams))
    names(res) <- nams
    
    i <- 0L
    info <- NULL # handles case where nams is empty, when no predictors
    ## For each predictor that "move" call Predict separately, and rbind results
    callenv <- parent.frame()
    for(nam in nams) {
      i <- i + 1L
      m$name <- nam
      lv     <- eval(m, callenv)
      j <- attr(lv, 'info')
      if(i == 1L) info <- j
      else {
        info$varying <- c(info$varying, j$varying)
        info$adjust  <- c(info$adjust,  j$adjust)
      }
      attr(lv, 'info') <- NULL
      lv$.predictor. <- nam
      res[[nam]] <- lv
    }
    lv <- do.call('rbind.data.frame', res)
    class(lv) <- c('Predict', 'data.frame')
    attr(lv, 'info') <- info
    return(lv)
  }
    
  f <- sum(assume != 9)	##limit list to main effects factors
  parms  <- at$parms
  label  <- at$label
  values <- at$values
  yunits <- fit$units
  units  <- at$units
  scale  <- fit$scale.pred
  if(! length(scale)) scale <- "X * Beta"

  if(! length(fun)) {
    ylab <- scale[1]
    if(length(time))
      ylab <- ylabPlotmath <- ylabhtml <-
        if(loglog) paste("log[-log S(", format(time), ")]", sep="")
        else paste(format(time), yunits, "Survival Probability")
    else if(scale[1] == 'X * Beta') {
      ylabPlotmath <- expression(X*hat(beta))
      ylabhtml     <- 'X&\beta;'
      }
    else ylabPlotmath <- ylabhtml <- ylab
  }
  else ylab <- ylabPlotmath <- ylabhtml <- ''
  
  if(ref.zero & length(time))
    stop("ref.zero=TRUE does not make sense with time given")
  
  if(fname == '') factors <- dotlist ## name not an argument
  else {
    factors <- list()
    for(g in fname) factors[[g]] <- NA
  }
  nf   <- length(factors)
  fnam <- names(factors)
  
  if(nf < 1) stop("must specify predictors to vary")
  
  which <- charmatch(fnam, name, 0L)
  if(any(which == 0L))
    stop(paste("predictors(s) not in model:",
               paste(names(factors)[which == 0L], collapse=" ")))
  
  if(any(assume[which] == 9L))
    stop("cannot plot interaction terms")

  lim <- Getlim(at, allow.null=TRUE, need.all=FALSE)
 
  fnam   <- names(factors)
  nf     <- length(factors)
  xadjdf <- lim$limits[2L, , drop=FALSE]
  xadj   <- unclass(xadjdf)
  varying <- NULL
  
  if(nf == 0L) return(as.data.frame(xadj))
  if(nf < f) { ## number of arguments < total number of predictors
    ## Some non-varying predictors
    settings <- xadj
    if(adj.zero) for(x in names(settings)) {
      i <- match(x, name)
      settings[[x]] <- if(assume[i] %in% c(5L, 8L)) parms[[i]][1]
      else if(length(V <- lim$values[[x]]) & is.character(V)) V[1]
      else 0L
    }
    for(n in fnam) settings[[n]] <- factors[[n]]
  }
  else settings <- factors
  
  for(i in 1L : nf) {
    n <- fnam[i]
    v <- settings[[n]]
    lv <- length(v)
    if(lv == 0L) stop('a predictor has zero length')
    if(lv == 1L && is.na(v))
      settings[[n]] <- value.chk(at, which(name == n), NA, np, lim)
    if(length(settings[[n]]) > 1L) varying <- c(varying, n)
  }
  if(prod(sapply(settings,length)) > 1e5)
    stop('it is not sensible to predict more than 100,000 combinations')
  settings <- expand.grid(settings)
  if(length(ioff)) settings[[offsetVariableName]] <- offset[[1]]
  adjust <- NULL
  for(n in name[assume != 9L & name %nin% fnam])
    adjust <- paste(adjust, n, "=", 
                    if(is.factor(xadj[[n]])) as.character(xadj[[n]])
                    else format(xadj[[n]]), " ", sep="")
  
  j <- assume != 9L
  label  <- label[j]
  units  <- units[j]
  assume <- assume[j]
  names(label) <- names(units) <- names(assume) <- name[j]
  at <- list(label=label, units=units, assume.code=assume)
  
  info <- list(varying=varying, adjust=adjust, Design=at,
               ylabPlotmath=ylabPlotmath, ylabhtml=ylabhtml,
               ylab=ylab, yunits=yunits,
               ref.zero=ref.zero, adj.zero=adj.zero, time=time,
               conf.int=conf.int)
  if(type == 'model.frame') {
    attr(settings, 'info') <- info
    return(settings)
  }
  ## Number of non-slopes
  nrp     <- num.intercepts(fit)
  nrpcoef <- num.intercepts(fit, 'coef')
  if(nrp > 0L && (kint < 1L || kint > nrp))
    stop('illegal intercept number for kint')

  beta <- fit$coefficients
  bootdone <- length(boot.Coef <- fit$boot.Coef) && usebootcoef
  if(bootdone && conf.type %in% c('individual','simultaneous')) {
    warning('bootstrap estimates ignored when conf.type="simultaneous" or "individual"')
    bootdone <- FALSE
  }
  isMean <- length(fun) && ! is.function(fun) && fun == 'mean'
  # Claude Sonnet 5 2026-11-20          9 lines
  # For a yde-fitted orm() object, fun='mean' is handled below via a
  # direct predict(type='mean', se.fit=) call, which gets a real
  # delta-method CI with no need for bootstrap reps or posterior draws
  # -- unlike the ordinary-fit mechanism these two guards exist for
  # (Mean.lrm()'s own conf.int code needs X, which isn't available at
  # the plain xb <- m(xb) call site below, hence bootcov/bayes being
  # required there). isormyde is exempted from both; fun=Mean(object)
  # (a function, not the string) is instead stopped on explicitly just
  # below, since the string form is the one implemented for yde.
  if(isMean && ! bootdone & conf.int > 0 & ! bayes & ! isormyde)
    stop('specifying fun="mean" with conf.int > 0 does not make sense when not running bootcov (with coef.reps=TRUE)')
  if(isMean && isorm && conf.int > 0 && ! isormyde)
    stop("fun='mean' not implemented for orm models when confidence intervals are requested")
  if(isMean && isormyde && bootdone && conf.int > 0)
    stop("fun='mean' is not yet implemented together with bootcov(coef.reps=TRUE) for an orm object using yde/ytrans")
  # Claude Sonnet 5 2026-11-22
  # borm() has no bootcov(coef.reps=TRUE) concept and bayes is always
  # TRUE for it, so the first two guards above never fire for a borm
  # fit regardless of isbormyde -- no analogous guard needed there;
  # predict(type='mean', cint=) already gets a real HPD interval with
  # no bootstrap/draws requirement, same as isormyde's delta-method SE.
  if((isormyde || isbormyde) && length(fun) && is.function(fun))
    stop("for an orm()/borm() object using yde/ytrans, specify fun='mean' instead of ",
         "fun=Mean(object) (or omit fun= and call predict(object, newdata, ",
         "type='mean') directly)")

  # Claude Sonnet 5 2026-11-25          entire block
  # fun='quantile', paralleling fun='mean' immediately above, but
  # scoped to isormyde only for now -- the real computation
  # (ormYdeFittedMean()'s new type='quantile' branch) is a delta
  # method built on infoMxop()/the information matrix, which has no
  # borm() counterpart (a Bayesian fit has posterior draws instead);
  # borm()'s own version would need a different, per-draw-then-HPD
  # mechanism (not yet built -- see Quantile.borm()'s own TODO).
  # Ordinary (non-yde) fits are not given this string shorthand either:
  # Quantile(fit) already has a working conf.int there via the existing
  # (non-yde) band-inversion mechanism, reachable today by passing an
  # actual function built from it as fun= (confirmed working directly,
  # unlike yde, where that path is structurally unavailable -- see
  # ormYdeQuantile()'s own doc comment).
  isQuantile <- length(fun) && ! is.function(fun) && fun == 'quantile'
  # Claude Sonnet 5 2026-11-26
  # isbormyde now has a real implementation too (bormPredictYde()'s own
  # type='quantile', a per-posterior-draw analog of ormYdeFittedMean()'s
  # delta method -- see borm.r), so the earlier stop() for isbormyde is
  # gone; only truly unsupported fit types (ordinary orm/lrm, blrm/borm
  # without yde) still stop here, matching isMean's own scope just above.
  if(isQuantile && ! isormyde && ! isbormyde)
    stop("fun='quantile' is currently only implemented for a yde-fitted ",
         "orm() or borm() object; for other fit types build fun from ",
         "Quantile(object) directly, e.g. ",
         "function(lp, ...) Quantile(object)(quantile, lp, ...)")
  if(isQuantile && bootdone && conf.int > 0)
    stop("fun='quantile' is not yet implemented together with ",
         "bootcov(coef.reps=TRUE) for an orm object using yde/ytrans")

  if(! length(time)) {
    xx <- Pred()
    if(length(attr(xx, "strata")) && any(is.na(attr(xx, "strata"))))
      warning("Computed stratum NA.  Requested stratum may not\nexist or reference values may be illegal strata combination\n")
    
    if(length(xx) == 0L)
      stop("model has no covariables and survival not plotted")
    xb <- if(is.list(xx)) xx$linear.predictors else xx
    if(isMean) {
      if(isormyde) {
        # Claude Sonnet 5 2026-11-20          entire block
        # Predict() already has settings (the actual covariate grid)
        # in scope, so it calls predict(type='mean') directly rather
        # than going through Mean()/ormYdeMean() -- no need to pass
        # settings out to a closure and have the closure call predict()
        # again. xx$lower/xx$upper are overwritten here (mean scale)
        # so the generic "lower <- as.vector(xx$lower)" step further
        # below picks these up unchanged; they otherwise still held
        # ormPredictYde()'s own type='lp' (log-odds-scale) CI from
        # Pred(), which would be the wrong scale for fun='mean'.
        r <- predict(fit, settings, type='mean', se.fit = conf.int > 0)
        if(conf.int > 0) {
          xb <- r$mean
          zcrit <- qnorm((1 + conf.int) / 2)
          xx$lower <- xb - zcrit * r$se.fit
          xx$upper <- xb + zcrit * r$se.fit
        } else xb <- r
      } else if(isbormyde) {
        # Claude Sonnet 5 2026-11-22
        # Same reasoning as isormyde just above, but bormPredictYde()'s
        # own type='mean' already returns an HPD interval directly
        # (list(mean=, lower=, upper=), via posterior draws) when
        # cint>0, rather than a delta-method se.fit Predict() would
        # then need to turn into bounds itself -- so lower/upper are
        # used as returned, no qnorm step needed here.
        r <- predict(fit, settings, type='mean',
                     posterior.summary=posterior.summary, cint=conf.int)
        if(conf.int > 0) {
          xb <- r$mean
          xx$lower <- r$lower
          xx$upper <- r$upper
        } else xb <- r
      } else if(isblrm || isborm) {
        # Claude Sonnet 5 2026-11-22
        # Nothing left to do: xb (and, via xx$lower/xx$upper, the CI)
        # were already transformed by fun=Mean(fit) inside Pred()'s own
        # call above (isormyde/isbormyde, which do need the separate
        # handling above, were already routed there). Calling m(xb)
        # again here would double-apply Mean() to an already-mean
        # value.
      } else {
        m <- Mean(fit)
        xb <- m(xb)
      }
      }
    # Claude Sonnet 5 2026-11-25/26          entire block
    # fun='quantile', structured exactly like the isMean block just
    # above (Predict() already has settings in scope, so it calls
    # predict(type='quantile') directly rather than going through
    # Quantile()/ormYdeQuantile()/bormYdeQuantile()). q=quantile
    # (Predict()'s own new formal, since a bare fun='quantile' string
    # cannot also carry which quantile is wanted the way fun=Mean(fit)
    # needs no such extra argument). xx$lower/xx$upper are overwritten
    # here (quantile scale) for the same reason isMean's own branch
    # overwrites them -- they otherwise still hold Pred()'s own
    # type='lp' (log-odds-scale) CI, the wrong scale for a quantile.
    # isormyde's delta-method SE and isbormyde's HPD interval are both
    # already on the quantile scale when predict() returns them (unlike
    # isMean's isormyde branch, which needs an extra zcrit*se.fit step
    # for the CI -- type='quantile' returns lower/upper directly from
    # both ormYdeFittedMean() and bormPredictYde(), so no such step is
    # needed here for either fit type).
    if(isQuantile) {
      if(isormyde) {
        r <- predict(fit, settings, type='quantile', q=quantile, cint=conf.int)
      } else if(isbormyde) {
        r <- predict(fit, settings, type='quantile', q=quantile,
                     posterior.summary=posterior.summary, cint=conf.int)
      }
      if(conf.int > 0) {
        xb <- r$quantile
        xx$lower <- r$lower
        xx$upper <- r$upper
      } else xb <- r
    }
    if(bootdone && conf.int > 0) {
      X <- Pred(type='x')
      pred <- t(matxv(X, boot.Coef,
                      kint=kint,  bmat=TRUE))
      if(isMean) {
        for(k in 1L : nrow(pred))
          pred[k,] <- m(pred[k,], intercepts=boot.Coef[k, 1L : nrp])
      }
      lim <- bootBCa(xb, pred, type=boot.type, n=nobs(fit), seed=fit$seed,
                     conf.int=conf.int)
      if(! is.matrix(lim)) lim <- as.matrix(lim)
      xx$lower <- lim[1L, ]
      xx$upper <- lim[2L, ]
    }    # end if(bootdone)
  }    # if(! length(time))
  else {   # time specified
    if(bootdone)
      stop('time may not be specified if bootcov was used with coef.reps=TRUE')
    xx <- survest(fit, settings, times=time, loglog=loglog, 
                  conf.int=conf.int)
    xb <- as.vector(xx$surv)
  }        # end time specified
  if(conf.int > 0) {
    lower <- as.vector(xx$lower)
    upper <- as.vector(xx$upper)
  }
  
  # Claude Sonnet 5 2026-11-17          2 lines
  # borm() also applies fun= itself (via predict.borm(), delegating to
  # predict.blrm() for the non-yde case) -- exclude it here too, same
  # reason blrm() already was, to avoid applying fun twice.
  if(! isblrm && ! isborm && length(fun) && is.function(fun)) {
    ## If fun is for example the result of Mean.lrm or Quantile.orm
    ## and confidence limits are requested, must use entire design matrix
    ## to get variances.  Note that conf.int must also have been requested
    ## when calling Mean/Quantile
    xb <- if(conf.int > 0 &&
             all(c('X', 'conf.int') %in% names(formals(fun)))) {
            X  <- Pred(type='x')
            fun(xb, X=X, conf.int=conf.int)
          } else fun(xb)
    if(conf.int > 0 && length(lims <- attr(xb, 'limits'))) {
      lower <- lims$lower
      upper <- lims$upper
      } else if(conf.int > 0) {
        lower <- fun(lower)
        upper <- fun(upper)
      }	
  }
  
  settings$yhat   <- xb
  if(conf.int > 0) {
    settings$lower  <- lower
    settings$upper  <- upper
  }
  class(settings) <- c('Predict', 'data.frame')
  attr(settings, 'info') <- info
  settings
}

print.Predict <- function(x, ...) {
  print.data.frame(x)
  info <- attr(x, 'info')
  cat('\nResponse variable (y):', info$ylab,'\n')
  if(length(info$adjust) == 1)
    cat('\nAdjust to:',info$adjust,'\n')
  ci <- info$conf.int
  if(ci > 0)
    cat('\nLimits are', ci, 'confidence limits\n')
  invisible()
}

perimeter <- function(x, y, xinc=diff(range(x))/10., n=10.,
                      lowess.=TRUE) {

  s <- ! is.na(x+y)
  x <- x[s]
  y <- y[s]
  m <- length(x)
  if(m < n) stop("number of non-NA x must be >= n")

  i <- order(x)
  x <- x[i]
  y <- y[i]
  s <- n : (m - n + 1L)
  x <- x[s]
  y <- y[s]

  x <- round(x / xinc) * xinc

  g <- function(y, n) {
    y <- sort(y)
    m <- length(y)
    if(n > (m - n + 1L)) c(NA, NA)
    else c(y[n], y[m - n + 1L])
  }
  
  r <- unlist(tapply(y, x, g, n=n))
  i <- seq(1L, length(r), by=2)
  rlo <- r[i]
  rhi <- r[-i]
  s <- ! is.na(rlo + rhi)
  if(! any(s))
    stop("no intervals had sufficient y observations")

  x <- sort(unique(x))[s]
  rlo <- rlo[s]
  rhi <- rhi[s]
  if(lowess.) {
    rlo <- lowess(x, rlo)$y
    rhi <- lowess(x, rhi)$y
  }

  structure(cbind(x, rlo, rhi),
            dimnames=list(NULL, c("x","ymin","ymax")),
            class='perimeter')
}

rbind.Predict <- function(..., rename) {
  d <- list(...)
  ns <- length(d)
  if(ns == 1) return(d[[1]])
  
  info <- attr(d[[1L]], 'info')
  
  if(! missing(rename)) {
    trans <- function(input, rename) {
      k <- input %in% names(rename)
      if(any(k)) input[k] <- rename[input[k]]
      input
    }
    info$varying <- trans(info$varying, rename)
    names(info$Design$label) <- trans(names(info$Design$label), rename)
    names(info$Design$units) <- trans(names(info$Design$units), rename)
    names(info$Design$assume.code) <-
      trans(names(info$Design$assume.code), rename)
  }
    
  info$Design$label <- c(info$Design$label, .set.='Set')
  info$Design$units <- c(info$Design$units, .set.='')
  info$varying <- c(info$varying, '.set.')

  sets <- names(d)
  if(! length(sets)) sets <- paste('Set', 1 : ns)
  obs.each.set <- sapply(d, function(x) length(x[[1]]))
  .set. <- rep(sets, obs.each.set)
  .set. <- factor(.set., levels=unique(.set.))

  info$adjust <- sapply(d, function(x) attr(x, 'info')$adjust)

  ## If first varying variable is not always the same but the second
  ## is, take varying[1] to be ".x."

  ## What in the heck is this for???

  if(FALSE) {
  first  <- sapply(d, function(x) attr(x, 'info')$varying[1])
  second <- sapply(d, function(x) {
    y <- attr(x, 'info')$varying
    if(length(y) < 2) '' else y[2] } )
  if((length(unique(first)) > 1) && (all(second == second[1])))
    info$varying[1] <- '.x.'
  }

  if(! missing(rename)) for(i in 1L : ns)
    names(d[[i]]) <- trans(names(d[[i]]), rename)

  result <- do.call('rbind.data.frame', d)
  result$.set. <- .set.
  attr(result, 'info') <- info
  class(result) <- c('Predict', 'data.frame')
  result
}
