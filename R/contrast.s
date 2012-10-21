contrast <- function(fit, ...) UseMethod("contrast")

contrast.rms <-
  function(fit, a, b, cnames=NULL,
           type=c('individual','average','joint'),
           conf.type=c('individual','simultaneous'),
           weights='equal', conf.int=0.95, tol=1e-7, expand=TRUE, ...)
{
  type <- match.arg(type)
  conf.type <- match.arg(conf.type)
  if(conf.type == 'simultaneous') require(multcomp)
  
  zcrit <- if(length(idf <- fit$df.residual)) qt((1+conf.int)/2, idf) else
  qnorm((1+conf.int)/2)

  da <- do.call('gendata', list(fit, factors=a, expand=expand))
  xa <- predict(fit, da, type='x')
  ma <- nrow(xa)

  if(missing(b))
    {
      xb <- 0*xa
      db <- da
    }
  else
    {
      db <- do.call('gendata', list(fit, factors=b, expand=expand))
      xb <- predict(fit, db, type='x')
    }
  mb <- nrow(xb)
  
  vary <- NULL
  if(type!='average' && !length(cnames))
    {
      ## If two lists have same length, label contrasts by any variable
      ## that has the same length and values in both lists
      if(ma==mb)
        {
          if(ncol(da) != ncol(db)) stop('program logic error')
          if(any(sort(names(da)) != sort(names(db))))
            stop('program logic error')
          k <- integer(0)
          nam <- names(da)
          for(j in 1:length(da))
            if(all(as.character(da[[nam[j]]])==as.character(db[[nam[j]]])))
              k <- c(k, j)
          if(length(k)) vary <- da[k]
        }
      else if(max(ma,mb)>1)
        {
          ## Label contrasts by values of longest variable in list if
          ## it has the same length as the expanded design matrix
          d <- if(ma>1) a else b
          l <- sapply(d, length)
          vary <- if(sum(l==max(ma,mb))==1)d[l==max(ma,mb)]
        }
    }
  
  if(max(ma,mb)>1 && min(ma,mb)==1)
    {
      if(ma==1) xa <- matrix(xa, nrow=mb, ncol=ncol(xb), byrow=TRUE)
      else
        xb <- matrix(xb, nrow=ma, ncol=ncol(xa), byrow=TRUE)
    }
  else if(mb != ma)
    stop('number of rows must be the same for observations generated\nby a and b unless one has one observation')

  X <- xa - xb
  p <- ncol(X)
  m <- nrow(X)
  
  if(is.character(weights))
    {
      if(weights!='equal') stop('weights must be "equal" or a numeric vector')
      weights <- rep(1, m)
    }
  else
    if(length(weights) > 1 && type != 'average')
      stop('can specify more than one weight only for type="average"')
    else
      if(length(weights) != m) stop(paste('there must be',m,'weights'))
  weights <- as.vector(weights)
  if(m > 1 && type=='average')
    X <- matrix(apply(weights*X, 2, sum) / sum(weights), nrow=1,
                dimnames=list(NULL,dimnames(X)[[2]]))
  
  est <- drop(X %*% coef(fit))
  v <- X %*% vcov(fit, regcoef.only=FALSE) %*% t(X)
  ndf <- if(is.matrix(v))nrow(v) else 1
  se <- if(ndf==1) sqrt(v) else sqrt(diag(v))
  Z <- est/se
  P <- if(length(idf)) 2*(1-pt(abs(Z), idf)) else 2*(1-pnorm(abs(Z)))
  if(conf.type=='individual') {
    lower <- est - zcrit*se
    upper <- est + zcrit*se
  } else {
    u <- confint(glht(fit, X,
                      df=if(length(idf)) idf else 0),
                      level=conf.int)$confint
    lower <- u[,'lwr']
    upper <- u[,'upr']
  }
  
  res <- list(Contrast=est, SE=se,
              Lower=lower, Upper=upper,
              Z=Z, Pvalue=P, 
              var=v, df.residual=idf,
              X=X, 
              cnames=if(type=='average')NULL else cnames,
              nvary=length(vary),
              conf.type=conf.type, conf.int=conf.int)
  if(type != 'average') res <- c(vary, res)
  
  r <- qr(v, tol=tol)
  nonred <- r$pivot[1:r$rank]   # non-redundant contrasts
  redundant <- (1:length(est)) %nin% nonred
  res$redundant <- redundant
  
  if(type=='joint')
    {
      est <- est[!redundant]
      v <- v[!redundant, !redundant, drop=FALSE]
      res$jointstat <- as.vector(est %*% solve(v, tol=tol) %*% est)
    }
  
  structure(res, class='contrast.rms')
}

print.contrast.rms <- function(x, X=FALSE, fun=function(u)u,
                               jointonly=FALSE, ...)
{
  edf <- x$df.residual
  sn <- if(length(edf))'t' else 'Z'
  pn <- if(length(edf))'Pr(>|t|)' else 'Pr(>|z|)'
  w <- x[1:(x$nvary + 6)]
  w$Z <- round(w$Z, 2)
  w$Pvalue <- round(w$Pvalue, 4)
  no <- names(w)
  no[no=='SE'] <- 'S.E.'
  no[no=='Z'] <- sn
  no[no=='Pvalue'] <- pn
  
  cnames <- x$cnames
  if(!length(cnames)) cnames <- if(x$nvary)rep('',length(x[[1]])) else
    as.character(1:length(x[[1]]))
  if(any(x$redundant)) cnames <- paste(ifelse(x$redundant, '*', ' '), cnames)
  attr(w,'row.names') <- cnames
  attr(w,'class') <- 'data.frame'
  w$Contrast <- fun(w$Contrast)
  w$SE       <- fun(w$SE)
  w$Lower    <- fun(w$Lower)
  w$Upper    <- fun(w$Upper)

  # Assign modified names to w
  names(w) <- no

  # Print w
  if(!jointonly)
    {
      print(as.matrix(w),quote=FALSE)
      if(any(x$redundant)) cat('\nRedundant contrasts are denoted by *\n')
    }
  
  jstat <- x$jointstat
  if(length(jstat))
    {
      cat('\nJoint test for all contrasts=0:\n\n')
      ndf <- sum(!x$redundant)
      if(length(edf))
        {
          Fstat <- jstat/ndf
          Pval <- 1 - pf(Fstat, ndf, edf)
          cat('F(',ndf,',',edf,')=',round(Fstat,3),', P=', round(Pval,4),
              '\n', sep='')
        }
      else
        {
          Pval <- 1 - pchisq(jstat, ndf)
          cat('Chi-square=', round(jstat,2),' with ', ndf, ' d.f.  P=',
              round(Pval,4),'\n', sep='')
        }
    }
  if(!jointonly && length(edf))cat('\nError d.f.=',edf,'\n')
  cat('\nConfidence intervals are', x$conf.int, x$conf.type,
      'intervals\n')
  if(X)
    {
      cat('\nDesign Matrix for Contrasts\n\n')
      if(is.matrix(x$X)) dimnames(x$X) <- list(cnames, dimnames(x$X)[[2]])
      print(x$X)
    }
  invisible()
}

