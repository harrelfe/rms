plot.xmean.ordinaly <- function(x, data, subset, na.action,
                                subn=TRUE, cr=FALSE, topcats=1,
                                cex.points=.75, ...)
{
  X <- match.call(expand.dots=FALSE)
  X$subn <- X$cr <- X$topcats <- X$cex.points <- X$... <- NULL
  if(missing(na.action)) X$na.action <- na.keep
  Terms <- if(missing(data)) terms(x)
  else
    terms(x, data=data)
  X$formula <- Terms
  X[[1]] <- as.name('model.frame')
  X <- eval.parent(X)
  resp <- attr(Terms, 'response')
  if(resp==0) stop('must have a response variable')
  
  nx <- ncol(X) - 1
  Y <- X[[resp]]
  nam <- as.character(attr(Terms, 'variables'))
  nam <- nam[-1]

  dopl <- function(x, y, cr, xname, yname)
    {
      s <- !is.na(unclass(Y)+x)
      y <- y[s]
      x <- x[s]
      n <- length(x)
      f <- lrm.fit(x, y)
      fy <- f$freq/n
      
      ##Following is pulled out of predict.lrm
      ns <- length(fy) - 1  # number of intercepts
      k <- ns + 1
      intcept <- f$coef[1:ns]
      xb <- f$linear.predictors - intcept[1]
      xb <- sapply(intcept, '+', xb)
      P <- 1/(1+exp(-xb))
      
      P <- matrix(P, ncol=ns)
      P <- cbind(1, P) - cbind(P, 0)  #one column per prob(Y=j)

      xmean.y <- tapply(x, y, mean)
      xp <- x*P/n
      xmean.y.po <- apply(xp, 2, sum)/fy
      yy <- 1:length(fy)
      rr <- c(xmean.y, xmean.y.po)
      if(cr)
        {
          u <- cr.setup(y)
          s <- u$subs
          yc <- u$y
          xc <- x[s]
          cohort <- u$cohort
          xcohort <- matrix(0, nrow=length(xc), ncol=length(levels(cohort))-1)
          xcohort[col(xcohort)==unclass(cohort)-1] <- 1  # generate dummies
          cof <- lrm.fit(cbind(xcohort, xc), yc)$coefficients
          cumprob <- rep(1, n)
          for(j in 1:k)
            {
              P[,j] <- cumprob* (if(j==k) 1
              else
               plogis(cof[1] + (if(j>1) cof[j] else 0) + cof[k]*x))
              cumprob <- cumprob - P[,j]
            }
          xp <- x*P/n
          xmean.y.cr <- apply(xp, 2, sum)/fy
          rr <- c(rr, xmean.y.cr)
        }
      plot(yy, xmean.y, type='b', ylim=range(rr),
           axes=FALSE, xlab=yname, ylab=xname, ...)
      mgp.axis(1, at=yy, labels=names(fy))
      mgp.axis(2)
      lines(yy, xmean.y.po, lty=2, ...)
      if(cr) points(yy, xmean.y.cr, pch='C', cex=cex.points)
      if(subn) title(sub=paste('n=',n,sep=''),adj=0)
    }
  

  for(i in 1:nx)
    {
      x <- X[[resp+i]]
      if(is.factor(x))
        {
          f <- table(x)
          ncat <- length(f)
          if(ncat < 2)
            {
              warning(paste('predictor',
                            nam[resp+i],'only has one level and is ignored'))
              next
            }
          nc <- min(ncat-1, topcats)
          cats <- (names(f)[order(-f)])[1:nc]
          for(wcat in cats)
            {
              xx <- 1*(x==wcat)
              xname <- paste(nam[resp+i], wcat, sep='=')
              dopl(xx, Y, cr, xname, nam[resp])
            }
        } else dopl(x, Y, cr, nam[resp+i], nam[resp])
    }
  invisible()
}
