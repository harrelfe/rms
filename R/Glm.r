#' rms Version of glm
#' 
#' This function saves `rms` attributes with the fit object so that
#' `anova.rms`, `Predict`, etc. can be used just as with `ols`
#' and other fits.  No `validate` or `calibrate` methods exist for
#' `Glm` though.
#' 
#' For the `print` method, format of output is controlled by the user
#' previously running `options(prType="lang")` where `lang` is
#' `"plain"` (the default), `"latex"`, or `"html"`.
#' 
#' 
#' @aliases Glm
#' @param formula,family,data,weights,subset,na.action,start,offset,control,model,method,x,y,contrasts
#' see [stats::glm()]; for `print` `x` is the result of `Glm`
#' @param ... ignored
#' @return a fit object like that produced by [stats::glm()] but with
#' `rms` attributes and a `class` of `"rms"`, `"Glm"`,
#' `"glm"`, and `"lm"`.  The `g` element of the fit object is
#' the \eqn{g}-index.
#' @seealso [stats::glm()],[Hmisc::GiniMd()], [prModFit()], [stats::residuals.glm]
#' @keywords models regression
#' @md
#' @examples
#' 
#' ## Dobson (1990) Page 93: Randomized Controlled Trial :
#' counts <- c(18,17,15,20,10,20,25,13,12)
#' outcome <- gl(3,1,9)
#' treatment <- gl(3,3)
#' f <- glm(counts ~ outcome + treatment, family=poisson())
#' f
#' anova(f)
#' summary(f)
#' f <- Glm(counts ~ outcome + treatment, family=poisson())
#' # could have had rcs( ) etc. if there were continuous predictors
#' f
#' anova(f)
#' summary(f, outcome=c('1','2','3'), treatment=c('1','2','3'))
#' 

Glm <- 
  function(formula, family = gaussian, data = environment(formula),
           weights, subset, na.action = na.delete, start = NULL, offset = NULL,
           control = glm.control(...), model = TRUE, method = "glm.fit",
           x = FALSE, y = TRUE, contrasts = NULL, ...)
{
  call <- match.call()
  if (is.character(family)) family <- get(family)
  if (is.function(family))  family <- family()
  if (!length(family$family)) {
    print(family)
    stop("`family' not recognized")
  }

  mt <- terms(formula, dta=data)

  callenv <- parent.frame()   # don't delay these evaluations
  weights <- if(! missing(weights)) eval(substitute(weights), data, callenv)
  subset  <- if(! missing(subset )) eval(substitute(subset),  data, callenv)

  mf <-
    modelData(data, formula,
              subset  = subset, weights=weights,
              na.action=na.action, callenv=callenv)

  mf <- Design(mf, formula=formula)

  at <- attributes(mf)
  desatr <- at$Design
  attr(mf, 'Design') <- NULL
  nact <- attr(mf, 'na.action')

  sformula   <- at$sformula
  mmcolnames <- desatr$mmcolnames
    
  switch(method, model.frame = return(mf), glm.fit = 1, 
         stop(paste("invalid `method':", method)))
  
  xvars <- as.character(attr(mt, "variables"))[-1]
  if ((yvar <- attr(mt, "response")) > 0)
    xvars <- xvars[-yvar]
    xlev <- if (length(xvars) > 0) {
      xlev <- lapply(mf[xvars], levels)
      xlev[!sapply(xlev, is.null)]
    }
  X <- if(! is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    intcpt <- if(attr(mt, 'intercept') > 0) '(Intercept)'
  alt <- attr(mmcolnames, 'alt')
  if(! all(mmcolnames %in% colnames(X)) && length(alt))
    mmcolnames <- alt
    X <- X[, c(intcpt, mmcolnames), drop=FALSE]
  colnames(X) <- c(if(length(intcpt)) 'Intercept', desatr$colnames)
  
#    colnames(X) <- if(attr(mt, 'intercept') > 0)
#    c('Intercept', desatr$colnames)
#    else desatr$colnames

  Y <- model.response(mf, "numeric")
  weights <- model.weights(mf)
  offset <- attr(mf, 'offset')
  if(!length(offset)) offset <- 0
  if (length(weights) && any(weights < 0))
    stop("Negative wts not allowed")
  if (length(offset) > 1 && length(offset) != NROW(Y))
    stop(paste("Number of offsets is", length(offset), ", should equal",
               NROW(Y), "(number of observations)"))
  fit <- glm.fit(x = X, y = Y, weights = weights, start = start,
                 offset = offset, family = family, control = control,
                 intercept = attr(mt, "intercept") > 0)
  if (length(offset) && attr(mt, "intercept") > 0) {
    fit$null.deviance <-
      if(is.empty.model(mt))
        fit$deviance
      else glm.fit(x = X[, "Intercept", drop = FALSE], y = Y,
                   weights = weights, start = start, offset = offset,
                   family = family, control = control,
                   intercept = TRUE)$deviance
  }
  if (model) fit$model <- mf
  if (x)  fit$x <- X[, -1, drop=FALSE]
  if (!y) fit$y <- NULL
  fit <- c(fit, list(call = call, formula = formula, sformula=sformula,
                     terms = mt, data = data, offset = offset,
                     control = control, method = method,
                     contrasts = attr(X, "contrasts"), xlevels = xlev,
                     Design=desatr, na.action=nact,
                     assign=DesignAssign(desatr,1,mt),
                     g=GiniMd(fit$linear.predictors)))
  class(fit) <- c('Glm', 'rms', 'glm', 'lm')
  fit
}
##' Print a `Glm` Object
##'
##' Prints a `Glm` object, optionally in LaTeX or html
##' @title print.glm
##' @param x `Glm` object
##' @param digits number of significant digits to print
##' @param coefs specify `coefs=FALSE` to suppress printing the table of
##' model coefficients, standard errors, etc.  Specify `coefs=n` to print
##' only the first `n` regression coefficients in the model.
##' @param title a character string title to be passed to `prModFit`
##' @param ... ignored
##' @author Frank Harrell
print.Glm <- function(x, digits=4, coefs=TRUE,
                      title='General Linear Model', ...)
{
  k <- 0
  z <- list()

  if(length(zz <- x$na.action))
    {
      k <- k + 1
      z[[k]] <- list(type=paste('naprint', class(zz)[1], sep='.'), list(zz))
    }

  cof <- coef(x)
  lr <- x$null.deviance - x$deviance
  dof <- x$rank - (names(cof)[1]=='Intercept')
  pval <- 1 - pchisq(lr, dof)

  ci <- x$clusterInfo
  misc <- reListclean(Obs=length(x$residuals),
                   'Residual d.f.'=x$df.residual,
                   'Cluster on'=ci$name,
                   Clusters=ci$n,
                   g = x$g)
  lr   <- reListclean('LR chi2'     = lr,
                   'd.f.'        = dof,
                   'Pr(> chi2)' = pval)
  headings <- c('', 'Model Likelihood\nRatio Test')
  data <-  list(c(misc, c(NA,NA,NA,NA,3)),
                c(lr,   c(2, NA,-4)))
  k <- k + 1
  z[[k]] <- list(type='stats', list(headings=headings, data=data))
  
  se <- sqrt(diag(vcov(x)))
  k <- k + 1
  z[[k]] <- list(type='coefmatrix',
                 list(coef=cof, se=se))
  prModFit(x, title=title, z, digits=digits, coefs=coefs, ...)
}

summary.Glm <- function(...) summary.rms(...)

vcov.Glm <- function(object, regcoef.only=TRUE, intercepts='all', ...) {
  v <- object$var
  if(!length(v)) v <- getS3method('vcov', 'glm')(object, ...)
  ns <- num.intercepts(object, 'var')
  if(ns > 0 && length(intercepts)==1 && intercepts=='none')
    v <- v[-(1 : ns), -(1 : ns), drop=FALSE]
  v
}

# Varcov.glm <- function(object, ...)
#{
#  if(length(object$var))
#    return(object$var)  ## for Glm
#  
#  s <- summary.glm(object)
#  s$cov.unscaled * s$dispersion
#}

##' Residuals for `Glm`
##'
##' This function mainly passes through to `residuals.glm` but for `type='score'` computes the matrix of score residuals using code modified from `sandwich::estfun.glm`.
##' @title residuals.Glm
##' @param object  a fit object produced by `Glm`
##' @param type either `'score'` or a `type` accepted by `residuals.glm`
##' @param ... ignored
##' @return a vector or matrix
##' @author Frank Harrell
residuals.Glm <- function(object, type, ...) {
  if(type == 'score') {
    if(! length(object[['x']])) stop('did not specify x=TRUE in fit')
    X <- cbind(Intercept=1, object$x)
    # Code modified from sandwich::estfun.glm
    w <- object$weights
    r <- object$residuals * w   # "working residuals"
    dispersion <- if(substring(object$family$family, 1, 17) %in%
      c('poisson', 'binomial', 'Negative Binomial')) 1 else sum(r ^2, na.rm=TRUE) / sum(w, na.rm=TRUE)
    return(r * X / dispersion)
  }
  residuals.glm(object, type=type, ...)
}

predict.Glm <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms","ccterms","adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE,
           conf.type=c('mean','individual','simultaneous'),
           kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               kint, na.action, expand.na, center.terms, ...)
  }

latex.Glm <- function(...) latexrms(...)
