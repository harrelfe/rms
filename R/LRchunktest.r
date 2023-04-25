LRchunktest <- function(object, i) {
  k   <- class(object)[1]
  fmi <- k == 'fit.mult.impute'
  if(fmi) k <- class(object)[2]
  
  X   <- object[['x']]
  y   <- object[['y']]
  if(! length(X) || ! length(y))
    stop('you must specify x=TRUE, y=TRUE in the fit to enable LR tests')
  
  if(length(i) == ncol(X)) {   # overall test used Model L.R.
    chisq <- object$stats['Model L.R.'] # and dealt with mult imputations
    return(c(chisq=chisq, df=ncol(X)))
  }
  
  X <- X[, - i, drop=FALSE]
  
  devf <- object$deviance
  devf <- devf[length(devf)]
  logl <- object$loglik
  logl <- logl[length(logl)]
  
  switch(k,
         lrm = {
            dev  <- lrm.fit(X, y, maxit=12, tol=1e-7)$deviance[2]
         },
         orm = {
           dev <- orm.fit(X, y, family=object$family)$deviance[2]
         },
         cph = {
           devf <- -2 * logl
           dev  <- -2 * coxphFit(X, y, method=object$method,
                            strata=object$strata,
                            type=attr(y, 'type'))$loglik[2]
         },
         psm = {
           devf <- -2 * logl
           dev  <- -2 * survreg.fit2(X, y, dist=object$dist,
                                     tol=1e-12)$loglik[2]
         },
         Glm = {
           dev <- glm.fit(x=cbind(1., X), y=y,
                          control=glm.control(),
                          offset=object$offset,
                          family=object$family)$deviance
         },
         stop('for LR tests, fit must be from lrm, orm, cph, psm, Glm')
         )
  chisq <- dev - devf
  if(fmi && length(object$fmimethod) && object$fmimethod != 'ordinary')
    chisq <- chisq / object$n.impute
  c(chisq=chisq, df=length(i))
  }
