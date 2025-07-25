# Quickly fits a model like an existing fit object but possibly changing
# X or y.  This is for models in rms using MLE.  For lrm and orm
# compstats is turned off by default.
# This function facilitates fast anova tables of LR tests,
# fast profile likelihood confidence intervals, and fast bootstrapping.
#
# For orm models, specify ytarget=constant to subset intercepts and reform
# the fit to be like a binary regression model for Y >= ytarget
# Set ytarget to NA to use the default reference intercept in doing this.
# The original fit (object) provides the default reference intercept and
# y values.
# Set storevals=FALSE to not store original X, y, offset, strata, weights as
# default argument values in generated functions
# storevals can be FALSE when the model re-fitting is specifying all data-related
# arguments to the fitter each time.  This applies when what='fitter'.
# If object is an orm and ytarget is used, set compvar=TRUE to compute the variance matrix
# that uses only the reference intercept.
# When storevals=TRUE, the subset argument in the generated function will be useful at times.
# When initial is provided to quickRefit(..., what='fitter') every fit run will use this initial value
# This is useful in bootstrapping and is used by orm.fit.  initial is used only for orm.fit.

quickRefit <-
  function(object,
           X              = object[['x']],
           y              = object[['y']],
           offset         = object[['offset'        ]],
           strata         = object[['strata'        ]],
           penalty.matrix = object[['penalty.matrix']],
           weights        = object[[if(k == 'Glm')'oweights' else 'weights']],
           compstats      = FALSE,
           gradtol        = 0.001,
           ytarget        = NULL,
           what           = c('chisq', 'deviance', 'fit', 'fitter'),
           compvar        = FALSE,
           initial        = NULL,
           storevals      = TRUE,
           ...) {

  what <- match.arg(what)
  k    <- class(object)[1]
  fmi  <- k == 'fit.mult.impute'
  if(fmi) k <- class(object)[2]
  if(k == 'ols' && length(penalty.matrix)) k <- 'olsp'
  if(k == 'orm' && length(ytarget)) k <- 'ormt'

  # is.null instead of ! length because X may be a no-covariate matrix during profile likelihood
  if(is.null(X) || is.null(y))    # is.null instead of ! length because X may be a no-covariate matrix during profile likelihood
    stop('you must specify x=TRUE, y=TRUE when fitting, or provide them to quickRefit')

  if(length(offset) && any(offset != 0e0) && (k %in% c('ols', 'olsp', 'Rq')))
    stop('offset not implemented for ols or Rq')

  # Need arguments strata, ytarget, penalty.matrix on all fitters so that predab.resample,
  # validate, calibrate, bootcov can pass them without a problem, even to fitters
  # that ignore them
  
g <-
  switch(k,
         ols = function()
                                   lm.fit.qr.bare(as.matrix(X)[subset,,drop=FALSE], y[subset], intercept=TRUE, ...),
         olsp= function()
                                   lm.pfit(cbind(Intercept=1e0, as.matrix(X)[subset,,drop=FALSE]), y[subset],
                                           penalty.matrix=penalty.matrix, regcoef.only=TRUE, ...),
         lrm = function()
           lrm.fit(as.matrix(X)[subset,,drop=FALSE], y[subset], compstats=compstats, offset=offset[subset], 
                   penalty.matrix=penalty.matrix, weights=weights[subset], opt_method=opt_method, ...),
         orm = function()
           orm.fit(as.matrix(X)[subset,,drop=FALSE], if(is.matrix(y)) y[subset,,drop=FALSE] else y[subset],
                        family=family, compstats=compstats, gradtol=gradtol, offset=offset[subset], initial=initial,
                        penalty.matrix=penalty.matrix, weights=weights[subset], opt_method=opt_method), #, ...),
         ormt= function() {
           f <- orm.fit(as.matrix(X)[subset,,drop=FALSE], if(is.matrix(y)) y[subset,,drop=FALSE] else y[subset],
                        family=family, compstats=compstats, gradtol=gradtol, offset=offset[subset], initial=initial,
                        penalty.matrix=penalty.matrix, weights=weights[subset], opt_method=opt_method) # $, ...)
           ns  <- f$non.slopes
           cof <- f$coefficients
           # For character y, find the intercept that exactly matchines ytarget
           # For numeric y, find the intercept corresponding to being closest to ytarget
           yu   <- f$yunique[-1]
           iref <- if(is.character(yu)) which(yu == ytarget) else which.min(abs(yu - ytarget))
           if(! length(iref)) stop('no intercept matches ytarget=', ytarget)
           i   <- c(iref, (ns + 1) : length(cof))
           cof <- cof[i]
           names(cof[1]) <- 'Intercept'
           attr(cof, 'intercepts') <- 1L
           f$coefficients          <- cof
           f$non.slopes            <- 1L
           if(compvar) f$var       <- infoMxop(f$info.matrix, i=i)
           f$info.matrix           <- NULL
           f
         } ,
         cph = function()
                                   coxphFit(as.matrix(X)[subset,,drop=FALSE], y[subset,,drop=FALSE], method=method, type=type,
                                            strata=strata[subset], weights=weights[subset], offset=offset[subset], ...),
         psm = function()
                                   survreg.fit2(as.matrix(X)[subset,,drop=FALSE], y[subset], offset=offset[subset], dist=dist,
                                                fixed=fixed, parms=parms, ...),
         Glm = function() {
           f <- glm.fit(x=cbind(1e0, as.matrix(X)[subset,,drop=FALSE]), y=as.vector(y)[subset], family=family, 
                        offset=offset[subset], weights=weights[subset],
                        control=control, ...)
           f$oweights <- weights
           f
         },
         Rq  = function()
                        quantreg::rq.wfit(cbind(Intercept=1e0, as.matrix(X)[subset,,drop=FALSE]), y[subset], tau=tau,
                                          weights=weights[subset], method=method, ...),
         bj  = function()
                    bj.fit(as.matrix(X)[subset,,drop=FALSE], y[subset], ...),
         stop('fit must be from ols, lrm, orm, cph, psm, Glm, bj, Rq')
         )

fm <- if(storevals) list(X=X,    y=y,    strata=strata, offset=offset, weights=weights, penalty.matrix=penalty.matrix)
      else          list(X=NULL, y=NULL, strata=NULL,   offset=NULL  , weights=NULL,    penalty.matrix=penalty.matrix)
formals(g) <- c(fm, alist(subset=TRUE, ytarget=, opt_method=, ...=))
#formals(g) <- c(formals(g), alist(subset=TRUE, ytarget=, opt_method=, ...))

  if(FALSE) {
  formals(g)$X       <- X
  formals(g)$y       <- y
  # Could not set arguments to NULL, as formals(g)$x <- NULL removes x as an argument
  formals(g)$strata  <- strata
  formals(g)$offset  <- offset
  formals(g)$weights <- weights
  formals(g)$penalty.matrix <- penalty.matrix
  }

if(k %in% c('lrm', 'orm', 'ormt')) {
                           formals(g)$compstats <- compstats
                           formals(g)$opt_method <- 'NR'
                           if(k != 'lrm') formals(g)$gradtol <- gradtol
}
if(k %in% c('orm', 'ormt')) {
                           fm <- list(initial=initial)
                           formals(g) <- c(formals(g), fm)
}
if(k == 'ormt' && is.na(ytarget)) {
  yu   <- object$yunique[-1]
  ytarget <- if(is.character(yu)) yu[object$interceptRef] else median(object[['y']])
}
if(k == 'ormt')            formals(g)$ytarget   <- ytarget
if(k == 'ormt')            formals(g)$compvar   <- compvar                                                    
if(k == 'cph')             formals(g)$type      <- attr(y, 'type')
if(k %in% c('cph', 'Rq'))  formals(g)$method    <- object$method
if(k == 'Rq')              formals(g)$tau       <- object$tau
if(k == 'psm') {
  formals(g)$dist           <- object$dist
  fixed <- object$fixed
  fixed <- if(length(fixed) == 1 && is.logical(fixed) && ! fixed) list()
            else list(scale=TRUE)
  formals(g)$fixed <- fixed
  formals(g)$parms <- object[['parms']]
  }
if(k %in% c('orm', 'ormt', 'Glm')) formals(g)$family <- object$family
if(k == 'Glm') formals(g)$control <- glm.control()
if(getOption('rmsdebug', FALSE)) {
  len <- function(x) if(is.matrix(x)) paste0(nrow(x), 'x', ncol(x))
           else if(length(x) == 1) as.character(x) else length(x)
  cat('\nArguments of function constructed by quickRefit:\n\n')
  print(sapply(formals(g), len), quote=FALSE)
}
if(what == 'fitter') return(g)
f <- g(X, y, ...)
if(what == 'fit' || (length(f$fail) && f$fail)) return(f)

dev  <- getDeviance(f, k)
fdev <- dev[length(dev)]
if(what == 'deviance') return(fdev)

  dev0 <- dev[1]
  dev  <- dev[length(dev)]
  if(what == 'deviance') return(dev)
  chisq <- dev0 - dev
  if(fmi && length(object$fmimethod) && object$fmimethod != 'ordinary')
    chisq <- chisq / object$n.impute
  c(chisq=chisq, df=ncol(X))
  }

# Specify fitclass manually if fit was from a quick shortcut fitter
# Only the first 3 letters are used from fitclass
# Or fitclass can be the whole original model fit object or all
# the classes from it
getDeviance <- function(object, fitclass) {
  if(missing(fitclass)) k <- class(object)
  else k <- if(is.character(fitclass)) fitclass else class(fitclass)
  if(k[1] == 'fit.mult.impute') k <- k[2] else k <- k[1]
  k <- substring(k, 1, 3)

  dname <- if(k %in% c('lrm', 'orm', 'Glm')) 'deviance' else
      if(k %in% c('cph', 'psm')) 'loglik' else 'notmle'
  if(dname == 'notmle') stop('fit did not use maximum likelihood estimation')

  dev <- object[[dname]]
  if(dname == 'loglik') dev <- -2e0 * dev
  if(k == 'Glm') dev <- c(object$null.deviance, dev)
  dev
}

utils::globalVariables(c('opt_method', 'method', 'dist', 'parms', 'control', 'tau'))
