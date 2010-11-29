hazard.ratio.plot <-
  function(x, Srv, which,
           times, e=30, subset, conf.int=.95, legendloc=NULL,
           smooth=TRUE, pr=FALSE, pl=TRUE, add=FALSE, ylim,cex=.5,
           xlab="t",ylab, antilog=FALSE, ...)
{
  if(missing(ylab)) ylab <- if(antilog)"Hazard Ratio" else "Log Hazard Ratio"

  trans <- if(antilog) function(x) exp(x) else function(x) x

  if(is.matrix(x))
    {
      nam <- dimnames(x)[[2]]
      if(!length(nam)) nam <- paste("x[",1:ncol(x),"]",sep="")
    }
  else
	{
      nam <- label(x)
      x <- as.matrix(oldUnclass(x))
      if(!length(nam)) nam <- ""
    }

  y <- Srv[,1];   event <- Srv[,2]
  if(length(y)!=nrow(x))stop("number of rows in x must be length of y")
  nx <- ncol(x)
  if(missing(which)) which <- 1:nx

  labele <- attr(Srv, "event.label")
  if(!length(labele)) labele <- ""

  isna <- is.na(matxv(x,rep(1,nx)) + y + event)

  if(!missing(subset))isna <- isna | (!subset)
  x <- x[!isna,,drop=FALSE]
  if(length(dimnames(x)[[2]])==0)
    dimnames(x) <- list(NULL,paste("x",1:nx,sep=""))
  y <- y[!isna]
  event <- event[!isna]

  if(!missing(times))uft<-c(0,sort(times),1000000)
  else
    {
      nblock<-max(round(sum(event)/e),2)
      uft<-c(0,quantile(y[event==1],
                        seq(0,1,length=nblock+1))[2:nblock], 1000000)
      uft <- unique(uft)
    }
  thr <- NULL
  lhr <- NULL
  se <- NULL
  for(i in seq(length(uft)-1))
    {
      s<-y>=uft[i]
      tt<-pmin(y[s],uft[i+1])
      ev<-event[s] & (y[s]<=uft[i+1])
      if(sum(ev)>nx)
        {
          cox <- coxphFit(x[s,,drop=FALSE], cbind(tt,ev),
                          iter.max=10, eps=.0001, method="efron",
                          type=attr(Srv, 'type'))
          if(!is.character(cox))
            {
              if(pr)
                {
                  r <- range(tt)
                  cat(paste("Time interval:",format(r[1]),"-",
                            format(r[2]),"  At Risk:",sum(s),
                            "  Events:",sum(ev),"\n"))
                  k <- cbind(cox$coefficients,sqrt(diag(cox$var)))
                  dimnames(k) <- list(names(cox$coefficients),
                                      c("Coef","S.E."))
                  print(k)
                }
              tmid <- mean(y[y>=uft[i] & y<=uft[i+1]])
              thr <- c(thr,tmid)
              lhr <- cbind(lhr,cox$coef)
              se <- cbind(se,sqrt(diag(cox$var)))
            }
        }
    }

  if(!pl) return(list(time=thr,log.hazard.ratio=lhr,se=se))

  zcrit<-qnorm((1+conf.int)/2)
  for(j in which)
    {
      lhrj <- lhr[j,]
      sej <- se[j,]
      labelx <- nam[j]
      if(missing(ylim)) ylim <- trans(range(c(lhrj+zcrit*sej,lhrj-zcrit*sej)))
      if(!add)
        {
          oldpar <- par(c('err','mar'))
          on.exit(par(oldpar))
          oldmar <- oldpar$mar
          if(labelx!="" & labele!="")oldmar[1]<-oldmar[1]+1
          par(err=-1,mar=oldmar)
          plot(thr,trans(lhrj),xlab=xlab,ylim=ylim,ylab=ylab,...)
        }
      else
        points(thr,trans(lhrj))
      lines(thr,trans(lhrj))
      lines(thr,trans(lhrj+zcrit*sej),lty=2)
      lines(thr,trans(lhrj-zcrit*sej),lty=2)
      leg <- c("Subset Estimate",paste(format(conf.int),"C.L."))
      ltype <- 1:2
      if(smooth & length(thr)>3)
        {
          lines(supsmu(thr, trans(lhrj)), lty=3)
          leg <- c(leg,"Smoothed")
          ltype <- c(ltype,3)
        }

      if(!add)
        {
          labels <- ""
          if(labelx != "")labels <- paste("Predictor:",labelx,"\n",sep="")
          if(labele != "")labels <- paste(labels,"Event:",labele,sep="")
          title(sub=labels,adj=1,cex=cex)

          if(!interactive() && !length(legendloc))legendloc <- "ll"
          if(!length(legendloc))
            {
              cat("Click left mouse button at upper left corner for legend\n")
              z <- locator(1)
              legendloc <- "l"
            }
          else
            if(legendloc[1]!="none")
              {
                if(legendloc[1]=="ll")
                  z <- list(x=par("usr")[1],y=par("usr")[3])
                else
                  z <- list(x=legendloc[1],y=legendloc[2]) 		  }	
          if(legendloc[1]!="none")legend(z,leg,lty=ltype,cex=cex,bty="n")
        }
    }
  list(time=thr,log.hazard.ratio=lhr,se=se)
}

