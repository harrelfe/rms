# Fast backward elimination using a slow but numerically stable version
# of the Lawless-Singhal method (Biometrics 1978), used in the SAS
# PHGLM and LOGIST procedures
# Uses function solvet, a slightly edited version of solve that passes
# the tol argument to qr.
# Modified 12Oct92 - if scale parameter present, ignore last row and col of cov
# Modified 22Sep93 - new storage format for design attributes
# Modified 1Mar94 - add k.aic
# Modified 4Mar96 - use S commands instead of avia if not under UNIX
# Modified 19Feb11 - added force argument
#
# F. Harrell 18Jan91

fastbw <- function(fit, rule=c("aic", "p"),
                   type=c("residual","individual","total"), sls=.05, aics=0, 
                   eps=1e-9, k.aic=2, force=NULL)
{
  rule <- match.arg(rule)
  type <- match.arg(type)
  
  ns <- num.intercepts(fit)
  if(length(force)) force <- force + ns
  L <- if(ns==0) NULL else 1:ns
  
  pt  <- length(fit$coef)
  p   <- pt - ns
  atr <- fit$Design
  
  assume <- atr$assume.code
  if(!length(assume)) stop("fit does not have design information")
  assign <- fit$assign
  nama <- names(assign)[1]
  asso <- 1*(nama=="(Intercept)" | nama=="Intercept")
  
  f <- sum(assume != 8)
  strt <- integer(f)
  len <- strt
  j <- 0
  for(i in 1:length(assume)) {
    if(assume[i] != 8) {
      j <- j+1
      aj <- assign[[j + asso]]
      strt[j] <- min(aj)
      len[j]  <- length(aj)
    }
  }
  name <- atr$name[assume != 8]
  ed <- as.integer(strt + len - 1)

  if(type == 'total') type <- 'residual'
  if(length(force) && type != 'individual')
    warning('force probably does not work unless type="individual"')

  factors.in <- 1:f
  parms.in   <- 1:pt

  ## Not needed if using solve() instead of avia
  ## Allocate work areas for avia
  ## s1 <- double(pt)
  ## s2 <- s1
  ## s3 <- double(2*pt)
  ## s4 <- s1
  ## vsub <- double(pt*pt)
  ## pivot <- integer(pt)

  factors.del <- integer(f)
  chisq.del   <- double(f)
  df.del      <- integer(f)
  resid.del   <- double(f)
  df.resid    <- integer(f)
  beta        <- fit$coef
  Cov  <- vcov(fit, regcoef.only=TRUE, intercepts='all')
  ## Above ignores scale parameters; 'all' for orm fits
  cov  <- Cov
  Coef <- matrix(NA, nrow=f, ncol=pt, dimnames=list(NULL, names(beta)))
  d    <- 0

  dor2 <- inherits(fit, 'ols') &&
    (length(fit$y) || (length(fit$fitted.values) &&
                       length(fit$residuals)))
  if(dor2) {
    ## X <- fit$x
    Y   <- if(length(fit$y))fit$y else fit$fitted.values + fit$residuals
    r2  <- double(f)
    sst <- sum((Y-mean(Y))^2)
    sigma2 <- fit$stats['Sigma']^2
    ## Get X'Y using b=(X'X)^-1 X'Y, X'X^-1 = var matrix / sigma2
    xpy <- matrix(solvet(Cov, beta, tol=eps)*sigma2, ncol=1)
    ypy <- sum(Y^2)
  }

  for(i in 1:f) {
    fi <- length(factors.in)
    ln <- len[factors.in]
    st <- as.integer(ns + c(1, 1 + cumsum(ln[-fi]))[1 : fi])
    en <- as.integer(st + ln - 1)
    if(any(en > nrow(cov))) stop('program logic error')
    crit.min       <- 1e10
    chisq.crit.min <- 1e10
    jmin  <- 0
    dfmin <- 0
    k     <- 0
    factors.in.loop <- factors.in  #indirect reference prob in S 3.1
    for(j in factors.in.loop) {
      k <- k + 1
      ## can't get this to work in R - CHECK:
      ##	z <- if(.R.)
      ##      .Fortran("avia",beta,cov,chisq=double(1),length(beta),
      ##               st[k]:en[k],
      ##               ln[k],df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE,
      ##               PACKAGE="Design") else
      ##      .Fortran("avia",beta,cov,chisq=double(1),length(beta),
      ##               st[k]:en[k],
      ##               ln[k],df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE)
      ##	chisq <- z$chisq
      ##	df <- z$df
      
      ##replace previous 5 statements with following 3 to use slow method
      q <- st[k] : en[k]
      chisq <- if(any(q %in% force)) Inf else
               beta[q] %*% solvet(cov[q,q], beta[q], tol=eps)
      df <- length(q)
      
      crit <- switch(rule, aic=chisq-k.aic * df, p=pchisq(chisq, df))
      if(crit < crit.min) {
        jmin     <- j
        crit.min <- crit
        chisq.crit.min <- chisq
        df.min   <- df
      }	
    }
    
    factors.in <- factors.in[factors.in != jmin]
    parms.in <- parms.in[parms.in < strt[jmin] | parms.in > ed[jmin]]
    if(length(parms.in)==0) q <- 1:pt else q <- (1:pt)[-parms.in]
    
    ## if(under.unix && !.R.) {
    ## z <- if(.R.)
    ##  .Fortran("avia",fit$coef,Cov,chisq=double(1),
    ##           pt,q,as.integer(pt-length(parms.in)),
    ##           df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE,
    ##           PACKAGE="Design") else
    ##  .Fortran("avia",fit$coef,Cov,chisq=double(1),
    ##           pt,q,as.integer(pt-length(parms.in)),
    ##           df=integer(1),eps,vsub,s1,s2,s3,s4,pivot,NAOK=TRUE)
    ## resid <- z$chisq
    ## resid.df <- z$df
    ##}
    
    ##replace previous 5 statements with following 2 to use slow method
    resid <- fit$coef[q] %*% solvet(Cov[q,q], fit$coef[q], tol=eps)
    resid.df <- length(q)

    del <- switch(type,
           residual   = switch(rule, aic=resid - k.aic*resid.df <= aics,
                                     p=1 - pchisq(resid,resid.df) > sls),
           individual = switch(rule, aic = crit.min <= aics,
                                     p   = 1 - crit.min > sls)	)
    if(del) {
      d              <- d + 1
      factors.del[d] <- jmin
      chisq.del  [d] <- chisq.crit.min
      df.del     [d] <- df.min
      resid.del  [d] <- resid
      df.resid   [d] <- resid.df
      if(length(parms.in)) {
        cov.rm.inv <- solvet(Cov[-parms.in, -parms.in], tol=eps)
        cov.cross  <- Cov[parms.in, -parms.in, drop=FALSE]
        w    <- cov.cross %*% cov.rm.inv
        beta <- fit$coef[parms.in] - w %*% fit$coef[-parms.in]
        cov  <- Cov[parms.in, parms.in] - w %*% t(cov.cross)
        cof  <- rep(0, pt)
        cof[parms.in] <- beta
        Coef[d,]      <- cof
        if(dor2) {
          ## yhat <- matxv(X[,parms.in,drop=F], beta)
          ## r2[d] <- 1 - sum((yhat-Y)^2)/sst
          ## sse = Y'(I - H)Y, where H = X*inv(X'X)*X'
          ##     = Y'Y - Y'X*inv(X'X)*X'Y
          ##     = Y'Y - Y'Xb
          sse <- ypy - t(xpy[parms.in, , drop=FALSE])%*%beta
          r2[d] <- 1 - sse/sst
        }
      }
      else {
        beta <- NULL; cov <- NULL
        if(dor2) r2[d] <- 0
      }
    }
    else break
  }
  
  if(d > 0) {
    u <- 1:d
    fd <- factors.del[u]
    if(dor2) {
      r2 <- r2[u]
      Coef <- Coef[u,, drop=FALSE]
    }
    res <- cbind(chisq.del[u], df.del[u],
                 1 - pchisq(chisq.del[u], df.del[u]),
                 resid.del[u], df.resid[u],
                 1 - pchisq(resid.del[u], df.resid[u]),
                 resid.del[u] - k.aic *   df.resid[u])
    labs <- c("Chi-Sq", "d.f.", "P", "Residual", "d.f.", "P", "AIC")
    dimnames(res) <- list(name[fd], labs)
    if(length(fd)==f) fk <- NULL else fk <- (1:f)[-fd]
  }
  else {
    fd <- NULL
    res <- NULL
    fk <- 1:f
  }
  
  nf <- name[fk]
  
  pd <- NULL
  if(d > 0) for(i in 1:d) pd <- c(pd, (strt[fd[i]] : ed[fd[i]]))

  if(length(fd) == f) fk <- NULL
  else
    if(d==0) fk <- 1:f
    else fk <- (1:f)[-fd]
  if(length(pd)==p) pk <- L
  else
    if(d==0) pk <- 1:pt
    else
      pk <- (1:pt)[-pd]
  
  if(length(pd) != p) {
    beta <- as.vector(beta)
    names(beta) <- names(fit$coef)[pk]
    dimnames(cov) <- list(names(beta),names(beta))
  }
  
  if(dor2) res <- cbind(res, R2=r2)
  r <- list(result=res, names.kept=nf, factors.kept=fk,
            factors.deleted=fd,
            parms.kept=pk, parms.deleted=pd, coefficients=beta, var=cov,
            Coefficients=Coef,
            force=if(length(force)) names(fit$coef)[force])
  class(r) <- "fastbw"
  r
}

		
print.fastbw <- function(x, digits=4, estimates=TRUE,...)
{

  res <- x$result
  fd <- x$factors.deleted
  if(length(fd)) {
    cres <- cbind(dimnames(res)[[1]], format(round(res[,1], 2)),
                  format(res[,2]),
                  format(round(res[,3], 4)), format(round(res[,4], 2)),
                  format(res[,5]), format(round(res[,6], 4)),
                  format(round(res[,7], 2)),
                  if(ncol(res) > 7)format(round(res[,8], 3)))
    dimnames(cres) <- list(rep("", nrow(cres)),
                           c("Deleted", dimnames(res)[[2]]))
    cat("\n")
    if(length(x$force))
      cat('Parameters forced into all models:\n',
          paste(x$force, collapse=', '), '\n\n')
    print(cres, quote=FALSE)
    if(estimates && length(x$coef)) {
      cat("\nApproximate Estimates after Deleting Factors\n\n")
      cof <- coef(x)
      vv <- if(length(cof)>1) diag(x$var) else x$var
      z <- cof/sqrt(vv)
      stats <- cbind(cof, sqrt(vv), z, 1 - pchisq(z^2,1))
      dimnames(stats) <- list(names(cof), c("Coef","S.E.","Wald Z","P"))
      print(stats, digits=digits)
    }
  }
  else cat("\nNo Factors Deleted\n")
  cat("\nFactors in Final Model\n\n")
  nk <- x$names.kept
  if(length(nk))print(nk, quote=FALSE)
  else cat("None\n")
}
