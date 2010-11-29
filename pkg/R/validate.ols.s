validate <-
  function(fit,  method="boot", B=40,
           bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
           pr=FALSE,...)
  UseMethod("validate")


#Resampling optimism of discrimination and reliability of an ols regression
#B: # reps
#bw=T to incorporate backward stepdown (using fastbw) with params rule,type,sls
#pr=T to print results of each bootstrap rep
#Requires: predab.resample, fastbw, ols
#Frank Harrell 11 June 91

validate.ols <- function(fit, method="boot",
	B=40,bw=FALSE,rule="aic",type="residual",
	sls=.05,aics=0,pr=FALSE, u=NULL, rel=">", tolerance=1e-7, ...)
{
  fit.orig <- fit
  
  penalty.matrix <- fit.orig$penalty.matrix
  
  discrim <- function(x, y, fit, iter, evalfit=FALSE, u=NULL, rel=NULL,
                      pr=FALSE, ...)
	{
      resid <- if(evalfit) fit$residuals else y - x

      n <- length(resid)
      sst <- (n-1)*var(y)   # sum(y^2) - (sum(y)^2)/n
      mse <- sum(resid^2)
      rsquare <- 1 - mse/sst
      mse <- mse/n

      if(evalfit)
        {	#Fit being examined on sample used to fit
          intercept <- 0
          slope     <- 1
        }
      else
        {
          if(length(fit$coef)==1) {intercept <- mean(y)-mean(x); slope <- 1}
          else
            {
              coef <- lsfit(x, y)$coef   #Note x is really x*beta from other fit
              intercept <- coef[1]
              slope     <- coef[2]
            }
        }

      z <- c(rsquare, mse, GiniMd(slope*x), intercept, slope)
      nam <- c("R-square", "MSE", "g", "Intercept", "Slope")
      if(length(u))
        {
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
  
  ols.fit <- function(x,y,tolerance=1e-7,backward, 
                      penalty.matrix=NULL, xcol=NULL, ...)
    {
      if(!length(x))
        {
          ybar <- mean(y)
          n <- length(y)
          residuals <- y - ybar
          v <- sum(residuals^2)/(n-1)
          return(list(coef=ybar, var=v/n, residuals=residuals, fail=FALSE))
        }
      if(length(penalty.matrix) > 0)
        {
          if(length(xcol))
            {
              xcol <- xcol[-1]-1   # remove position for intercept
              penalty.matrix <- penalty.matrix[xcol,xcol,drop=FALSE]
            }
          fit <- lm.pfit(x, y, penalty.matrix=penalty.matrix,
                         tol=tolerance)
        }
      else
        {
          fit <- lm.fit.qr.bare(x,as.vector(y),tolerance=tolerance,
                                intercept=FALSE, xpxi=TRUE)
          if(backward) 
            fit$var <- sum(fit$residuals^2)*fit$xpxi/
              (length(y) - length(fit$coefficients))
        }
      c(fit,fail=FALSE)
    }
  
  predab.resample(fit.orig,method=method,fit=ols.fit,measure=discrim,pr=pr,
                  B=B,bw=bw,rule=rule,type=type,sls=sls,aics=aics,
                  tolerance=tolerance,
                  backward=bw,u=u, penalty.matrix=penalty.matrix,
                  rel=rel, ...)
}

print.validate <- function(x, digits=4, ...)
  print(round(unclass(x), digits), ...)
