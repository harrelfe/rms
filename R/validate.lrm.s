#Resampling optimism of discrimination and reliability of a logistic 
#regression model
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep

validate.lrm <- function(fit,method="boot",
	B=40, bw=FALSE, rule="aic", type="residual",
	sls=.05, aics=0, force=NULL, estimates=TRUE, pr=FALSE,
    kint,
	Dxy.method=if(k==1)"somers2" else "lrm",
	emax.lim=c(0,1), ...)
{
  k <- fit$non.slopes
  y <- fit$y
  if(length(y)==0) stop("fit did not use x=TRUE,y=TRUE")
  if(!is.factor(y)) y <- factor(y)   ## was category 11Apr02
  fit$y <- oldUnclass(y)-1  #mainly for Brier score (B)
  
  if(missing(kint)) kint <- floor((k+1)/2)
  
  penalty.matrix <- fit$penalty.matrix
  
  discrim <- function(x, y, fit, iter, evalfit=FALSE, pr=FALSE,
                      Dxy.method="somers2",
                      penalty.matrix, kint, ...)
    {
      if(evalfit) {	# Fit was for bootstrap sample
        stats <- fit$stats
        lr <- stats["Model L.R."]
        Dxy <- if(Dxy.method=="lrm") stats["Dxy"] else
        somers2(x,y)["Dxy"]
        intercept <- 0
        shrink <- 1
        n  <- stats["Obs"]
        D  <- (lr - 1)/n
        U  <- -2/n
        Q  <- D - U
        R2 <- stats["R2"]
        g  <- stats['g']
        gp <- stats['gp']
      }
      else {	
        k <- fit$non.slopes
        null.model <- length(fit$coefficients)==k
        refit <- if(null.model) lrm.fit(y=y) else lrm.fit(x,y,tol=1e-13)
        kr <- refit$non.slopes
        ## Model with no variables = null model
        stats <- refit$stats
        lr <- stats["Model L.R."]
        Dxy <- if(Dxy.method=="lrm") stats["Dxy"] else
        somers2(x,y)["Dxy"]
        intercept <- refit$coefficients[kint]
        shrink <- if(null.model) 1 else refit$coefficients[kr + 1]
        n <- stats["Obs"]
        D <- (lr-1)/n
        L01 <- -2 * sum( (y >= kint)*x - logb(1 + exp(x)), na.rm=TRUE)
        U <- (L01 - refit$deviance[2] - 2)/n
        Q <- D - U
        R2 <- stats["R2"]
        g  <- GiniMd(shrink*x)
        gp <- GiniMd(plogis(intercept + shrink*x))
      }
      P <- plogis(x)  # 1/(1+exp(-x))
      B <- sum(((y >= kint) - P)^2)/n
      z <- c(Dxy, R2, intercept, shrink, D, U, Q, B, g, gp)
      names(z) <- c("Dxy", "R2", "Intercept", "Slope", "D", "U", "Q", "B",
                    "g",   "gp")
      z
    }
  
  lrmfit <- function(x, y, maxit=12, tol=1e-7, penalty.matrix=NULL, 
                     xcol=NULL, ...)
    {
      if(length(xcol) && length(penalty.matrix) > 0)
        penalty.matrix <- penalty.matrix[xcol, xcol, drop=FALSE]
      lrm.fit(x, y, maxit=maxit, penalty.matrix=penalty.matrix, tol=tol)
    }

  z <- predab.resample(fit, method=method, fit=lrmfit, measure=discrim, pr=pr,
                       B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                       force=force, estimates=estimates, Dxy.method=Dxy.method,
                       non.slopes.in.x=FALSE,
                       penalty.matrix=penalty.matrix, kint=kint, ...)
  kept <- attr(z, 'kept')
  calib <- z[3:4,5]
  p <- seq(emax.lim[1],emax.lim[2],.0005)
  L <- logb(p/(1-p))
  P <- plogis(calib[1]+calib[2]*L)  # 1/(1+exp(-calib[1]-calib[2]*L))
  emax <- max(abs(p-P), na.rm=TRUE)
  z <- rbind(z[1:4,],c(0,0,emax,emax,emax,z[1,6]),z[5:nrow(z),])
  dimnames(z) <- list(c("Dxy", "R2","Intercept", "Slope", "Emax", "D", "U", "Q",
                        "B", "g", "gp"),
                      c("index.orig","training","test","optimism",
                        "index.corrected","n"))
  structure(z, class='validate', kept=kept)
}


validate.orm <- function(fit, method="boot",
	B=40, bw=FALSE, rule="aic", type="residual",
	sls=.05, aics=0, force=NULL, estimates=TRUE, pr=FALSE,  ...)
{
  k <- fit$non.slopes
  y <- fit$y
  if(length(y)==0) stop("fit did not use x=TRUE, y=TRUE")
  if(!is.factor(y)) y <- factor(y)
  
  penalty.matrix <- fit$penalty.matrix
  
  discrim <- function(x, y, fit, iter, evalfit=FALSE, pr=FALSE,
                      penalty.matrix, ...)
    {
      if(evalfit) {	 # Fit was for bootstrap sample
        stats <- fit$stats
        lr  <- stats["Model L.R."]
        rho <- stats["rho"]
        shrink <- 1
        n   <- stats["Obs"]
        R2  <- stats["R2"]
        g   <- stats['g']
        pdm <- stats['pdm']
      }
      else {
        k <- fit$non.slopes
        null.model <- length(fit$coefficients)==k
        refit <- if(null.model) ormfit2(y=y) else ormfit2(x, y, tol=1e-13)
        kr <- refit$non.slopes
        ## Model with no variables = null model
        stats <- refit$stats
        lr <- stats["Model L.R."]
        rho <- stats['rho']
        shrink <- if(null.model) 1 else refit$coefficients[kr + 1]
        n  <- stats["Obs"]
        R2 <- stats["R2"]
        g  <- GiniMd(shrink*x)
        pdm <- stats['pdm']
      }
      z <- c(rho, R2, shrink, g, pdm)
      names(z) <- c("rho", "R2", "Slope", "g", "pdm")
      z
    }
  
  ormfit2 <- function(x, y, maxit=12, tol=1e-7, penalty.matrix=NULL, 
                      xcol=NULL, ...)
    {
      if(length(xcol) && length(penalty.matrix) > 0)
        penalty.matrix <- penalty.matrix[xcol, xcol, drop=FALSE]
      # x has names() like y >= ... - DROP ??
      # predab.resample is getting constant x as if an intercept
      orm.fit(x, y, maxit=maxit, penalty.matrix=penalty.matrix, tol=tol)
    }

  z <- predab.resample(fit, method=method, fit=ormfit2, measure=discrim, pr=pr,
                       B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                       force=force, estimates=estimates, 
                       non.slopes.in.x=FALSE,
                       penalty.matrix=penalty.matrix,
                       allow.varying.intercepts=TRUE, ...)
  kept <- attr(z, 'kept')
  dimnames(z) <- list(c("rho", "R2", "Slope", "g", "pdm"),
                      c("index.orig","training","test","optimism",
                        "index.corrected","n"))
  structure(z, class='validate', kept=kept)
}
