calibrate.default <- function(fit, predy, 
			      method=c("boot","crossvalidation",".632","randomization"),
			      B=40, bw=FALSE, rule=c("aic","p"),
			      type=c("residual","individual"),
			      sls=.05, aics=0, force=NULL, pr=FALSE, kint,
			      smoother="lowess", digits=NULL, ...)
{
  call   <- match.call()
  method <- match.arg(method)
  rule   <- match.arg(rule)
  type   <- match.arg(type)

  ns <- num.intercepts(fit)
  if(missing(kint)) kint <- floor((ns+1)/2)
  clas <- attr(fit,"class")
  model <- if(any(clas=="lrm"))"lr"
           else if(any(clas=="ols")) "ol"
           else stop("fit must be from lrm or ols")
  lev.name <- NULL
  yvar.name <- as.character(formula(fit))[2]
  y <- fit$y
  n <- length(y)
  if(length(y)==0) stop("fit did not use x=TRUE,y=TRUE")
  if(model=="lr")
    {
      y <- factor(y)
      lev.name <- levels(y)[kint+1]
      fit$y <- as.integer(y)-1
    }

  predicted <- if(model=="lr") 
                 1/(1+exp(-(fit$linear.predictors-fit$coefficients[1] +
                            fit$coefficients[kint])))
               else
                 fit$linear.predictors
 
  if(missing(predy))
    {
      if(n < 11) stop("must have n > 10 if do not specify predy")
      p <- sort(predicted)
      predy <- seq(p[5], p[n-4], length=50)
      p <- NULL
    }

  penalty.matrix <- fit$penalty.matrix

  cal.error <- function(x, y, iter, smoother, predy, kint, model,
                        digits=NULL, ...)
    {
      if(model=="lr")
        {
          x <- plogis(x)
          y <- y >= kint
        }
      if(length(digits)) x <- round(x, digits)
      smo <- if(is.function(smoother)) smoother(x, y) else
       lowess(x, y, iter=0)
      cal <- approx(smo, xout=predy, ties=function(x)x[1])$y
      if(iter==0) storeTemp(cal,".orig.cal")
      cal-predy
    }

  fitit <- function(x, y, model, penalty.matrix=NULL, xcol=NULL, ...)
    {
    if(length(penalty.matrix) && length(xcol))
      {
        if(model=='ol') xcol <- xcol[-1] - 1   # take off intercept position
        penalty.matrix <- penalty.matrix[xcol,xcol,drop=FALSE]
      }
    switch(model,
           lr=lrm.fit(x, y, penalty.matrix=penalty.matrix, tol=1e-13),
           ol=c(if(length(penalty.matrix)==0)
             {
                  w <- lm.fit.qr.bare(x, y, intercept=FALSE, xpxi=TRUE)
                  w$var <- w$xpxi * sum(w$residuals^2) /
                    (length(y) - length(w$coefficients))
                  w
                }
                else 
                  lm.pfit(x, y, penalty.matrix=penalty.matrix), fail=FALSE))
  }

  z <- predab.resample(fit, method=method, fit=fitit, measure=cal.error,
                       pr=pr, B=B, bw=bw, rule=rule, type=type, sls=sls,
                       aics=aics, force=force,
                       non.slopes.in.x=model=="ol",
                       smoother=smoother, predy=predy, model=model, kint=kint,
                       penalty.matrix=penalty.matrix, ...)

  z <- cbind(predy, calibrated.orig=.orig.cal,
             calibrated.corrected=.orig.cal - z[,"optimism"],
             z)
  structure(z, class="calibrate.default", call=call, kint=kint, model=model,
            lev.name=lev.name, yvar.name=yvar.name, n=n, freq=fit$freq,
            non.slopes=ns, B=B, method=method, 
            predicted=predicted, smoother=smoother)
}

print.calibrate.default <- function(x, B=Inf, ...)
{
  at <- attributes(x)
  cat("\nEstimates of Calibration Accuracy by ",at$method," (B=",at$B,")\n\n",
      sep="")
  dput(at$call)
  if(at$model=="lr")
    {
      lab <- paste("Pr{",at$yvar.name,sep="")
      if(at$non.slopes==1) lab <- paste(lab,"=",at$lev.name,"}",sep="")
      else lab <- paste(lab,">=",at$lev.name,"}",sep="")
    }
  else lab <- at$yvar.name
  
  cat("\nPrediction of",lab,"\n\n")
  predicted <- at$predicted
  if(length(predicted))
    {  ## for downward compatibility
      s <- !is.na(x[,'predy'] + x[,'calibrated.corrected'])
      err <- predicted - approx(x[s,'predy'],x[s,'calibrated.corrected'], 
                                xout=predicted, ties=mean)$y
      cat('\nn=',length(err),    '   Mean absolute error=',
          round(mean(abs(err),na.rm=TRUE),3),'   Mean squared error=',
          round(mean(err^2,na.rm=TRUE),5),
          '\n0.9 Quantile of absolute error=',
          round(quantile(abs(err),.9,na.rm=TRUE),3),	   '\n\n',sep='')
    }
  print.default(x)
  kept <- at$kept
  if(length(kept))
    {
      cat("\nFactors Retained in Backwards Elimination\n\n")
      varin <- ifelse(kept, '*', ' ')
      print(varin[1:min(nrow(varin), B),], quote=FALSE)
      cat("\nFrequencies of Numbers of Factors Retained\n\n")
      nkept <- apply(kept, 1, sum)
      tkept <- table(nkept)
      names(dimnames(tkept)) <- NULL
      print(tkept)
    }
  
  invisible()
}

plot.calibrate.default <- function(x, xlab, ylab, xlim, ylim, legend=TRUE, 
                                   subtitles=TRUE, scat1d.opts=NULL, ...)
{
  at <- attributes(x)
  if(missing(ylab))
    ylab <- if(at$model=="lr") "Actual Probability"
    else paste("Observed", at$yvar.name)
  
  if(missing(xlab))
    {
      if(at$model=="lr")
        {
          xlab <- paste("Predicted Pr{",at$yvar.name,sep="")
          if(at$non.slopes==1)
            {
              xlab <- if(at$lev.name=="TRUE") paste(xlab, "}", sep="")
              else paste(xlab,"=", at$lev.name, "}", sep="")
            }
          else xlab <- paste(xlab,">=", at$lev.name, "}", sep="")
        }
      else xlab <- paste("Predicted", at$yvar.name)
    }
  
  p     <- x[,"predy"]
  p.app <- x[,"calibrated.orig"]
  p.cal <- x[,"calibrated.corrected"]
  if(missing(xlim) & missing(ylim))
    xlim <- ylim <- range(c(p, p.app, p.cal), na.rm=TRUE)
  else
    {
      if(missing(xlim)) xlim <- range(p)
      if(missing(ylim)) ylim <- range(c(p.app, p.cal, na.rm=TRUE))
    }
  
  plot(p, p.app, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, type="n", ...)
  predicted <- at$predicted
  err <- NULL
  if(length(predicted))
    {  ## for downward compatibility
      s <- !is.na(p + p.cal)
      err <- predicted - approx(p[s], p.cal[s], xout=predicted, ties=mean)$y
      cat('\nn=',n <- length(err),    '   Mean absolute error=',
          round(mae <- mean(abs(err), na.rm=TRUE),3),'   Mean squared error=',
          round(mean(err^2, na.rm=TRUE),5),
          '\n0.9 Quantile of absolute error=',
          round(quantile(abs(err), .9, na.rm=TRUE),3),	   '\n\n', sep='')
      if(subtitles) title(sub=paste('Mean absolute error=', round(mae,3),
                            ' n=', n, sep=''), cex=.65, adj=1)
      do.call('scat1d', c(list(x=predicted), scat1d.opts))
    }
  
  lines(p, p.app, lty=3)
  lines(p, p.cal, lty=1)
  abline(a=0, b=1, lty=2)
  if(subtitles) title(sub=paste("B=", at$B, "repetitions,", at$method), adj=0)
  if(!(is.logical(legend) && !legend))
    {
      if(is.logical(legend)) legend <- list(x=xlim[1] + .55*diff(xlim),
                                            y=ylim[1] + .32*diff(ylim))
      legend(legend, c("Apparent", "Bias-corrected", "Ideal"),
             lty=c(3,1,2), bty="n")
    }
  invisible(err)
}
