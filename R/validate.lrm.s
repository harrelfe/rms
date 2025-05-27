#Resampling optimism of discrimination and reliability of a logistic
#regression model
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep

validate.lrm <-
  function(fit,method="boot",
	         B=40, bw=FALSE, rule="aic", type="residual",
	         sls=.05, aics=0, force=NULL, estimates=TRUE, pr=FALSE,
           kint,
	         Dxy.method,
	         emax.lim=c(0,1), ...)
{
  if(! missing(Dxy.method)) warning('Dxy.method is deprecated and ignored')
  k <- fit$non.slopes
  y <- fit$y
  if(length(y)==0) stop("fit did not use x=TRUE,y=TRUE")
  if(!is.factor(y)) y <- factor(y)   ## was category 11Apr02
  fit$y <- unclass(y) - 1  #mainly for Brier score (B)

  if(missing(kint)) kint <- floor((k+1)/2)

  penalty.matrix <- fit$penalty.matrix

  discrim <- function(x, y, fit, iter, evalfit=FALSE, pr=FALSE,
                      penalty.matrix, kint, ...)
    {
      k <- fit$non.slopes
      null.model <- length(fit$coefficients) == k
      if(evalfit) {	# Fit was for bootstrap sample
        stats     <- fit$stats
        lr        <- stats["Model L.R."]
        Dxy       <- stats["Dxy"]
        intercept <- if(null.model) NA else 0
        shrink    <- if(null.model) NA else 1
        Emax      <- 0
        n         <- stats["Obs"]
        D         <- (lr - 1)/n
        U         <- -2 / n
        Q         <- D - U
        R2        <- stats["R2"]
        g         <- stats['g']
        gp        <- stats['gp']
      }
      else {
        refit     <- if(null.model) lrm.fit(y=y) else lrm.fit(x, y)
        kr        <- refit$non.slopes
        ## Model with no variables = null model
        stats     <- refit$stats
        lr        <- stats["Model L.R."]
        Dxy       <- stats["Dxy"]
        intercept <- if(null.model) NA else refit$coefficients[kint]
        shrink    <- if(null.model) NA else refit$coefficients[kr + 1]
        Emax      <- if(null.model) NA else {
          p <- seq(emax.lim[1], emax.lim[2], 0.0005)
          L <- qlogis(p)
          P <- plogis(intercept + shrink * L)
          max(abs(p - P), na.rm=TRUE)
        }
        n         <- stats["Obs"]
        D         <- (lr - 1) / n
        xc        <- pmin(40e0, pmax(x, -40e0))
        L01       <- -2 * sum( (y >= kint) * xc - logb(1 + exp(xc)), na.rm=TRUE)
        U         <- (L01 - refit$deviance[2] - 2)/n
        Q         <- D - U
        R2        <- stats["R2"]
        g         <- if(null.model) 0 else GiniMd(shrink * x)
        gp        <- if(null.model) 0 else GiniMd(plogis(intercept + shrink * x))
      }
      P <- plogis(x)  # 1/(1+exp(-x))
      B <- sum(((y >= kint) - P)^2)/n
      z <- c(Dxy, R2, intercept, shrink, Emax, D, U, Q, B, g, gp)
      names(z) <- c("Dxy", "R2", "Intercept", "Slope", "Emax", "D", "U", "Q", "B",
                    "g",   "gp")
      z
    }

  lrmfit <- function(x, y, penalty.matrix=NULL,
                     xcol=NULL, strata, iter, ...)
    {
      if(length(xcol) && length(penalty.matrix) > 0)
        penalty.matrix <- penalty.matrix[xcol, xcol, drop=FALSE]
      lrm.fit(x, y, penalty.matrix=penalty.matrix, ...)
    }

  z <- predab.resample(fit, method=method, fit=lrmfit, measure=discrim, pr=pr,
                       B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                       force=force, estimates=estimates,
                       non.slopes.in.x=FALSE,
                       penalty.matrix=penalty.matrix, kint=kint, ...)
  kept <- attr(z, 'kept')
  cn <- c("index.orig","training","test","optimism",
          "index.corrected","Lower", "Upper", "n")
  if('Lower' %nin% colnames(z)) cn <- cn[-(6:7)]

  dimnames(z) <- list(c("Dxy", "R2","Intercept", "Slope", "Emax", "D", "U", "Q",
                        "B", "g", "gp"), cn)
  structure(z, class='validate', kept=kept)
}


validate.orm <- function(fit, method="boot",
	B=40, bw=FALSE, rule="aic", type="residual",
	sls=.05, aics=0, force=NULL, estimates=TRUE, pr=FALSE,  ...)
{
  k <- fit$non.slopes
  y <- fit[['y']]
  if(length(y)==0) stop("fit did not use x=TRUE, y=TRUE")
  cens <- NCOL(y) == 2

  db <- getOption('validate.debug', FALSE)

  discrim <- function(x, y, fit, iter, evalfit=FALSE, pr=FALSE, ...)
    {
      if(evalfit) {	 # Fit was for bootstrap sample
        stats <- fit$stats
        lr  <- stats["Model L.R."]
        rho <- stats["rho"]
        Dxy <- stats["Dxy"]
        shrink <- 1
        n   <- stats["Obs"]
        R2  <- stats["R2"]
        g   <- stats['g']
        pdm <- stats['pdm']
      }
      else {
        k          <- fit$non.slopes
        null.model <- length(fit$coefficients) == k
        refit      <- if(null.model) ormfit2(y=y) else ormfit2(x, y)
        kr         <- refit$non.slopes
        ## Model with no variables = null model
        stats      <- refit$stats
        lr         <- stats["Model L.R."]
        rho        <- stats['rho']
        Dxy        <- stats['Dxy']
        shrink     <- if(null.model) NA else refit$coefficients[kr + 1]
        n          <- stats["Obs"]
        R2         <- stats["R2"]
        g          <- if(null.model) 0 else GiniMd(shrink * x)
        pdm        <- stats['pdm']
      }
      z <- c(rho, Dxy, R2, shrink, g, pdm)
      names(z) <- c("rho", "Dxy", "R2", "Slope", "g", "pdm")
      if(cens) z <- z[-1]   # rho doesn't handle censoring, is NA
      z
    }

  ormfit2 <- quickRefit(fit, what='fitter', storevals=FALSE, compstats=TRUE)
  if(db) {
    saveRDS(ormfit2, '/tmp/ormfit2.rds')
    ormfit3 <- function(...) {
      saveRDS(list(...), '/tmp/ormfit3.rds')
      ormfit2(...)
    }
  }

  z <- predab.resample(fit, method=method, fit=if(db) ormfit3 else ormfit2, measure=discrim,
                       pr=pr, B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                       force=force, estimates=estimates,
                       non.slopes.in.x=FALSE,
                       allow.varying.intercepts=TRUE, ...)
  kept <- attr(z, 'kept')
  cn <- c("index.orig","training","test","optimism",
          "index.corrected","Lower", "Upper", "n")
  if('Lower' %nin% colnames(z)) cn <- cn[-(6:7)]
  dimnames(z) <- list(c(if(! cens) "rho", "Dxy", "R2", "Slope", "g", "pdm"), cn)
  structure(z, class='validate', kept=kept)
}
