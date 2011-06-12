ols <- function(formula, data, weights, subset, na.action=na.delete, 
                method = "qr", model = FALSE, x = FALSE, y = FALSE,
                se.fit=FALSE, linear.predictors=TRUE,
                penalty=0, penalty.matrix, tol=1e-7, sigma=NULL,
                var.penalty=c('simple','sandwich'), ...)
{
  call <- match.call()
  var.penalty <- match.arg(var.penalty)
  m <- match.call(expand = FALSE)
  mc <- match(c("formula", "data", "subset", "weights", "na.action"), 
              names(m), 0)
  m <- m[c(1, mc)]
  m$na.action <- na.action
  m$drop.unused.levels <- TRUE
  m[[1]] <- as.name("model.frame")
  ##X's present) 
  if(length(attr(terms(formula),"term.labels")))
    {
      ## R's model.frame.default gives wrong model frame if [.factor
      ## removes unused factor levels
      dul <- .Options$drop.unused.levels
      if(!length(dul) || dul)
        {
          on.exit(options(drop.unused.levels=dul))
          options(drop.unused.levels=FALSE)
        }

      X <- Design(eval.parent(m))
      options(drop.unused.levels=dul)
      atrx <- attributes(X)
      atr <- atrx$Design
      nact <- atrx$na.action
      Terms <- atrx$terms
      assig <- DesignAssign(atr, 1, Terms)
      
      penpres <- FALSE
      if(!missing(penalty) && any(unlist(penalty) != 0)) penpres <- TRUE
      if(!missing(penalty.matrix) && any(penalty.matrix != 0)) penpres <- TRUE
      
      if(penpres && missing(var.penalty))
        warning('default for var.penalty has changed to "simple"')
      
      if(method == "model.frame") return(X)
      scale <- as.character(formula[2])
      attr(Terms, "formula") <- formula
      weights <- model.extract(X, 'weights')
      if(length(weights) && penpres)
        stop('may not specify penalty with weights')

      Y <- model.extract(X, 'response')
      n <- length(Y)
      if(model) m <- X
      X <- model.matrix(Terms, X)
      if(length(atr$colnames)) 
        dimnames(X)[[2]] <- c("Intercept",atr$colnames)
      else dimnames(X)[[2]] <- c("Intercept",dimnames(X)[[2]][-1])
      if(method=="model.matrix") return(X)				   }
  
  ##Model with no covariables:
  
  else
    {
      if(length(weights))
        stop('weights not implemented when no covariables are present')
      assig <- NULL
      yy <- attr(terms(formula),"variables")[1]
      Y <- eval(yy,sys.parent(2))
      nmiss <- sum(is.na(Y))
      if(nmiss==0) nmiss <- NULL else names(nmiss) <- as.character(yy)
      Y <- Y[!is.na(Y)]
      yest <- mean(Y)
      coef <- yest
      n <- length(Y)
      if(!length(sigma)) sigma <- sqrt(sum((Y-yest)^2)/(n-1))
      cov <- matrix(sigma*sigma/n, nrow=1, ncol=1,
                    dimnames=list("Intercept","Intercept"))
      fit <- list(coefficients=coef, var=cov,
                  non.slopes=1, fail=FALSE, residuals=Y-yest,
                  df.residual=n-1, intercept=TRUE)
      if(linear.predictors)
        {
          fit$linear.predictors <- rep(yest,n); 
          names(fit$linear.predictors) <- names(Y)
        }
      if(model) fit$model <- m
      if(x) fit$x <- matrix(1, ncol=1, nrow=n, 
                            dimnames=list(NULL,"Intercept"))
      if(y) fit$y <- Y
      fit$fitFunction <- c('ols','lm')
      oldClass(fit) <- c("ols","rms","lm")
      return(fit)
    }
  
  if(!penpres)
    {
      fit <- if(length(weights))
        lm.wfit(X, Y, weights, method=method,  ...)
      else 
        lm.fit(X, Y, method=method, ...)
      cov.unscaled <- chol2inv(fit$qr$qr)
      r <- fit$residuals
      yhat <- Y - r
      if(length(weights))
        { ## see summary.lm
          sse <- sum(weights * r^2)
          m <- sum(weights * yhat / sum(weights))
          ssr <- sum(weights * (yhat - m)^2)
          r2 <- ssr / (ssr + sse)
          if(!length(sigma)) sigma <- sqrt(sse/fit$df.residual)
        }
      else
        {
          sse <- sum(fit$residuals^2)
          if(!length(sigma)) sigma <- sqrt(sse/fit$df.residual)
          r2 <- 1-sse/sum((Y-mean(Y))^2)
        }
      fit$var <- sigma*sigma*cov.unscaled
      cnam <- dimnames(X)[[2]]
      dimnames(fit$var) <- list(cnam, cnam)
      fit$stats <- c(n=n,'Model L.R.'=-n*logb(1-r2),
                     'd.f.'=length(fit$coef)-1, R2=r2, g=GiniMd(yhat),
                     Sigma=sigma)
    }
  else
    {
      p <- length(atr$colnames)
      if(missing(penalty.matrix)) penalty.matrix <- Penalty.matrix(atr,X)
      if(nrow(penalty.matrix)!=p || ncol(penalty.matrix)!=p) 
        stop('penalty matrix does not have',p,'rows and columns')
      psetup <- Penalty.setup(atr, penalty)
      penalty <- psetup$penalty
      multiplier <- psetup$multiplier
      if(length(multiplier)==1) penalty.matrix <- multiplier*penalty.matrix
      else
        {
          a <- diag(sqrt(multiplier))
          penalty.matrix <- a %*% penalty.matrix %*% a
        }
      fit <- lm.pfit(X, Y,
                     penalty.matrix=penalty.matrix, tol=tol,
                     var.penalty=var.penalty)
      fit$penalty <- penalty
    }
  
  if(model) fit$model <- m
  if(linear.predictors)
    {
      fit$linear.predictors <- Y-fit$residuals
      names(fit$linear.predictors) <- names(Y)
    }
  if(x) fit$x <- X
  if(y) fit$y <- Y
  if(se.fit)
    {
      se <- drop((((X %*% fit$var) * X) %*% rep(1, ncol(X)))^0.5)
      names(se) <- names(Y)
      fit$se.fit <- se
    }
  fit <- c(fit, list(call=call, terms=Terms, Design=atr,
                     non.slopes=1, na.action=nact,
                     scale.pred=scale, fail=FALSE,
                     fitFunction=c('ols','lm')))
  fit$assign <- assig
  oldClass(fit) <- c("ols","rms","lm")
  fit
}


lm.pfit <- function(X, Y, penalty.matrix, tol=1e-7, regcoef.only=FALSE,
                    var.penalty=c('simple','sandwich'))
{
  var.penalty <- match.arg(var.penalty)
  p <- ncol(X)-1
  pm <- rbind(matrix(0, ncol=p+1, nrow=1),
              cbind(matrix(0, ncol=1, nrow=p), penalty.matrix))
  xpx <- t(X) %*% X
  Z <- solvet(xpx + pm, tol=tol)
  coef <- Z %*% t(X) %*% Y
  if(regcoef.only) return(list(coefficients=coef))
  yhat <- drop(X %*% coef)
  res  <- Y - yhat
  n <- length(Y)
  sse <- sum(res^2)
  s2 <- drop( (sse + t(coef) %*% pm %*% coef) / n )
  var <- if(var.penalty=='simple') s2 * Z else s2 * Z %*% xpx %*% Z
  cnam <- dimnames(X)[[2]]
  dimnames(var) <- list(cnam, cnam)
  sst <- (n-1)*var(Y)
  lr <- n*(1+logb(sst/n))-n*logb(s2)-sse/s2
  s2.unpen <- sse/n
  dag <- diag((xpx / s2.unpen) %*% (s2 * Z))
  df <- sum(dag) - 1
  stats <- c(n=n, 'Model L.R.'=lr, 'd.f.'=df, R2=1-sse/sst,
             g=GiniMd(yhat), Sigma=sqrt(s2))
  
  list(coefficients=drop(coef), var=var, residuals=res, df.residual=n-1,
       penalty.matrix=penalty.matrix, 
       stats=stats, effective.df.diagonal=dag)
}


predict.ols <- 
  function(object, newdata,
           type=c("lp","x","data.frame","terms","cterms","ccterms","adjto",
             "adjto.data.frame", "model.frame"),
           se.fit=FALSE, conf.int=FALSE, conf.type=c('mean','individual'),
           incl.non.slopes, non.slopes, kint=1,
           na.action=na.keep, expand.na=TRUE, center.terms=type=="terms", ...)
  {
    type <- match.arg(type)
    predictrms(object, newdata, type, se.fit, conf.int, conf.type,
               incl.non.slopes, non.slopes, kint,
               na.action, expand.na, center.terms, ...)
  }

print.ols <- function(x, digits=4, long=FALSE, coefs=TRUE, latex=FALSE,
                      title="Linear Regression Model", ...)
{
  k <- 0
  z <- list()
  

  if(length(zz <- x$na.action))
    {
      k <- k + 1
      z[[k]] <- list(type=paste('naprint', class(zz)[1], sep='.'), list(zz))
    }

  stats <- x$stats

  pen <- length(x$penalty.matrix) > 0

  resid <- x$residuals

  n <- length(resid)
  p <- length(x$coef)-(names(x$coef)[1]=="Intercept")
  if(length(stats)==0) cat("n=", n,"   p=",p,"\n\n",sep="")
  ndf <- stats['d.f.']
  df <- c(ndf, n-ndf-1, ndf)
  r2 <- stats['R2']
  sigma <- stats['Sigma']
  rdf <- df[2]
  rsqa <- 1 - (1 - r2)*(n-1)/rdf
  lrchisq <- stats['Model L.R.']
  if(lst <- length(stats))
    {
      misc <- reVector(Obs=stats['n'],
                       sigma=sigma,
                       'd.f.'=df[2])
      lr   <- reVector('LR chi2'     = lrchisq,
                       'd.f.'        = ndf,
                       'Pr(> chi2)' = 1 - pchisq(lrchisq, ndf))
      disc <- reVector(R2=r2, 'R2 adj'=rsqa, g=stats['g'])
      headings <- list('',
                       c('Model Likelihood', 'Ratio Test'),
                       c('Discrimination', 'Indexes'))
      data <- list(c(misc, c(NA,digits,NA)), c(lr, c(2,NA,4)), c(disc,3))
      k <- k + 1
      z[[k]] <- list(type='stats', list(headings=headings, data=data))
    }
  if(rdf > 5)
    {
      if(length(dim(resid)) == 2)
        {
          rq <- apply(t(resid), 1, quantile)
          dimnames(rq) <- list(c("Min", "1Q", "Median", "3Q",
                                 "Max"), dimnames(resid)[[2]])
        }
	  else
        {
          rq <- quantile(resid)
          names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        }
      k <- k + 1
      z[[k]] <- list(type=if(latex)'latexVector' else 'print',
                     list(rq, digits=digits),
                     tex=latex, title='Residuals')
    }
  else
    if(rdf > 0)
      {
        k <- k + 1
        z[[k]] <- list(type=if(latex)'latexVector' else 'print',
                       list(resid, digits=digits),
                       tex=latex, title='Residuals')
      }
  
  if(nsingular <- df[3] - df[1])
    {
      k <- k + 1
      z[[k]] <- list(type='cat',
                     paste(nsingular, 'coefficients not defined because of singularities'))
    }
      
  k <- k + 1
  se <- sqrt(diag(x$var))
  z[[k]] <- list(type='coefmatrix',
                 list(coef    = x$coefficients,
                      se      = se,
                      errordf = rdf))
  
  if(!pen)
    {
      if(long && p > 0)
        {
          correl <- diag(1/se) %*% x$var %*% diag(1/se)
          dimnames(correl) <- dimnames(x$var)
          cat("\nCorrelation of Coefficients:\n")
          ll <- lower.tri(correl)
          correl[ll] <- format(round(correl[ll], digits), ...)
          correl[!ll] <- ""
          k <- k + 1
          z[[k]] <- list(type='print', 
                         list(correl[-1,  - (p+1), drop = FALSE],
                              quote=FALSE, digits = digits))
        }
    }
  prModFit(x, title=title, z, digits=digits,
           coefs=coefs, latex=latex, ...)
  
  invisible()
}
