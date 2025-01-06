# Quickly fits a model like an existing fit object but possibly changing
# X or y.  This is for models in rms using MLE.  For lrm and orm
# compstats is turned off by default.
# This function facilitates fast anova tables of LR tests and
# fast profile likelihood confidence intervals
# Should also streamline bootcov to use this  TODO
#
# For orm models, specify ytarget=constant to subset intercepts and reform
# the fit to be like a binary regression model for Y >= ytarget
# Set ytarget to NA to use the default reference intercept in doing this.
# Set storevals=FALSE to not store original X, y, offset, strata, weights as
# default argument values in generated functions
# storevals can be FALSE when the model re-fitting is specifying all data-related
# arguments to the fitter each time.  This applies when what='fitter'.

quickRefit <-
  function(object,
           X              = object[['x']],
           y              = object[['y']],
           offset         = object[['offset'        ]],
           strata         = object[['strata'        ]],
           penalty.matrix = object[['penalty.matrix']],
           weights        = object[[if(k == 'Glm')'oweights' else 'weights']],
           compstats      = FALSE, 
           ytarget        = NULL,
           what           = c('chisq', 'deviance', 'fit', 'fitter'),
           storevals      = TRUE,
           ...) {

  what <- match.arg(what)
  k    <- class(object)[1]
  fmi  <- k == 'fit.mult.impute'
  if(fmi) k <- class(object)[2]
  if(k == 'ols' && length(penalty.matrix)) k <- 'olsp'
  if(k == 'orm' && length(ytarget)) {
    k <- 'ormt'
    if(what == 'fitter') stop('what ="fitter" is not allowed with ytarget')
  }
  if(! length(X) || ! length(y))
    stop('you must specify x=TRUE, y=TRUE when fitting, or provide them to quickRefit')

  if(length(offset) && any(offset != 0e0) && (k %in% c('ols', 'olsp', 'Rq')))
    stop('offset not implemented for ols or Rq')

  # Need arguments strata, ytarget, penalty.matrix on all fitters so that predab.resample,
  # validate, calibrate, bootcov can pass them without a problem, even to fitters
  # that ignore them
  
  g <-
  switch(k,
         ols = function(X, y, strata=NULL, offset=NULL, weights=NULL, penalty.matrix, ytarget, ...)
                                   lm.fit.qr.bare(X, y, intercept=TRUE, ...),
         olsp= function(X, y,  strata=NULL, offset=NULL, weights=NULL, penalty.matrix, ytarget, ...)
                                   lm.pfit(cbind(Intercept=1e0, X), y,
                                           penalty.matrix=penalty.matrix, regcoef.only=TRUE, ...),
         lrm = function(X, y, compstats, offset=NULL, penalty.matrix, weights=NULL,
                        strata=NULL, ytarget, ...)
           lrm.fit(X, y, compstats=compstats, offset=offset, 
                   penalty.matrix=penalty.matrix, weights=weights, ...),
         orm = function(X, y, family, compstats, offset=NULL, penalty.matrix, weights=NULL,
                        strata=NULL, ytarget, ...)
           orm.fit(X, y, family=family, compstats=compstats, offset=offset,
                        penalty.matrix=penalty.matrix, weights=weights,...),
         ormt= function(X, y, family, compstats, offset=NULL, penalty.matrix, weights=NULL,
                        strata=NULL, ytarget, ...) {
           f <- orm.fit(X, y, family=family, compstats=compstats, offset=offset,
                        penalty.matrix=penalty.matrix, weights=weights,
                        strata, ytarget, ...)
           ns  <- f$non.slopes
           cof <- f$coefficients
           if(! is.na(ytarget)) {
             # Y values corresponding to intercepts
             yu <- f$yunique[-1]
             # Linearly interpolate to return an intercept aimed
             # at Y >= ytarget
             intcept <- approx(yu, cof[1:ns], xout=ytarget)$y
             intattr <- approx(yu, 1:ns, xout=ytarget)$y
           } else {
             k         <- f$interceptRef
             intattr   <- k
             intercept <- cof[k]
           }
          names(intcept) <- 'Intercept'
          cof <- c(intcept, cof[(ns + 1) : length(cof)])
          attr(cof, 'intercepts') <- intattr
          f$coefficients          <- cof
          return(f)
         } ,
         cph = function(X, y, method, type, strata=NULL, weights=NULL, offset=NULL,
                        ytarget, penalty.matrix, ...)
                                   coxphFit(X, y, method=method, type=type,
                                            strata=strata, weights=weights, offset=offset, ...),
         psm = function(X, y, offset=NULL, dist, fixed, parms, strata=NULL, weights=NULL,
                        penalty.matrix, ytarget, ...)
                                   survreg.fit2(X, y, offset=offset, dist=dist,
                                                fixed=fixed, parms=parms, ...),
         Glm = function(X, y, family, offset=NULL, weights=NULL, control=glm.control(),
                        strata=NULL, penalty.matrix, ytarget, ...) {
           f <- glm.fit(x=cbind(1e0, X), y=as.vector(y), family=family, 
                        offset=offset, weights=weights,
                        control=control, ...)
           f$oweights <- weights
           f
         },
         Rq  = function(X, y, tau, offset=NULL, weights=NULL, method, strata=NULL,
                        penalty.matrix, ytarget,...)
                        quantreg::rq.wfit(cbind(Intercept=1e0, X), y, tau=tau,
                                          weights=weights, method=method, ...),
         bj  = function(X, y,  strata=NULL, offset=NULL, penalty.matrix, ytarget, ...)
                    bj.fit(X, y, ...),
         stop('fit must be from ols, lrm, orm, cph, psm, Glm, bj, Rq')
         )

if(storevals) {
  formals(g)$X       <- X
  formals(g)$y       <- y
  # Could not set arguments to NULL, as formals(g)$x <- NULL removes x as an argument
  formals(g)$strata  <- strata
  formals(g)$offset  <- offset
  formals(g)$weights <- weights
}

formals(g)$penalty.matrix <- penalty.matrix

if(k %in% c('lrm', 'orm', 'ormt'))
                           formals(g)$compstats <- compstats
if(k == 'ormt')            formals(g)$ytarget   <- ytarget
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
if(k %in% c('orm', 'Glm')) formals(g)$family <- object$family

if(what == 'fitter') return(g)
f <- g(X, y, ...)
if(what == 'fit') return(f)

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
