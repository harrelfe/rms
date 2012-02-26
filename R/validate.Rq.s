validate.Rq <-
  function(fit, method="boot",
           B=40, bw=FALSE, rule="aic", type="residual",
           sls=.05, aics=0, force=NULL,
           pr=FALSE, u=NULL, rel=">", tolerance=1e-7, ...)
{
  Rqfit <- RqFit(fit, wallow=FALSE)
  rqfit <-
    if(bw) function(x, y, ...) {  # need covariance matrix
      if(length(colnames(x)) && colnames(x)[1]=='Intercept')
        x <- x[,-1]
      w <- Rq(y ~ x, tau=fit$tau, method=fit$method, se=fit$se, hs=fit$hs)
      w$fail <- FALSE
      w
    }
    else {
      function(...) {
        w <- Rqfit(...)
        w$fail <- FALSE
        w
      }
    }
    
  fit.orig <- fit
  fit.orig$fail <- FALSE
  
  discrim <- function(x, y, fit, iter, evalfit=FALSE, u=NULL, rel=NULL,
                      pr=FALSE, ...)
	{
      resid <- if(evalfit) fit$residuals else y - x
      mad <- mean(abs(resid))
      if(evalfit) {	#Fit being examined on sample used to fit
        intercept <- 0
        slope     <- 1
      }
      else {
        if(length(fit$coef)==1) {intercept <- median(y)-mean(x); slope <- 1}
        else {
          cof <- Rqfit(cbind(1,x), y)$coefficients
          ##Note x is really x*beta from other fit
          intercept <- cof[1]
          slope     <- cof[2]
        }
      }
      z <- c(mad, cor(x, y, method='spearman'),
             GiniMd(slope*x), intercept, slope)
      nam <- c("MAD", "rho", "g", "Intercept", "Slope")
      if(length(u)) {
        yy <- if(rel==">") ifelse(y>u,  1, 0)
        else if(rel==">=") ifelse(y>=u, 1, 0)
        else if(rel=="<")  ifelse(y<u,  1, 0)
        else ifelse(y <= u, 1, 0)
        z <- c(z, somers2(x,yy)["Dxy"])
        nam <- c(nam, paste("Dxy Y",rel,format(u),sep=""))
        if(rel==">"|rel==">=") P <- 1-pnorm((u-x)/sqrt(mse))
        else P <- pnorm((u-x)/sqrt(mse))
        P0 <- sum(yy)/n
        L <- -2*sum(yy*logb(P)+(1-yy)*logb(1-P))
        L0<- -2*sum(yy*logb(P0)+(1-yy)*logb(1-P0))
        R2 <- (1-exp(-(L0-L)/n))/(1-exp(-L0/n))
        z <- c(z, R2)
        nam <- c(nam, paste("R2 Y",rel,format(u),sep=""))
      }
      names(z) <- nam
      z
    }
  
  predab.resample(fit.orig, method=method, fit=rqfit,
                  measure=discrim, pr=pr,
                  B=B, bw=bw, rule=rule, type=type, sls=sls, aics=aics,
                  force=force, tolerance=tolerance,
                  backward=bw, u=u, rel=rel, ...)
}
