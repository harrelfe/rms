#Function to divide x (e.g. Cox predicted survical at time y created by
#survest) into g quantile groups, get Kaplan-Meier estimates at time u
#(a scaler), and to return a matrix with columns x=mean x in
# quantile, n=#subjects, events=#events, and KM=K-M survival at time u,
# std.err = s.e. of log-log K-M
#Failure time=y   Censoring indicator=e
#Instead of supplying g, the user can supply the number of subjects to have
#in the quantile group on the average; then g will be computed.  The default
#m is 50, so the default g is n/50.
#If cuts is given (e.g. cuts=c(0,.1,.2,...,.9,.1)), it overrides m and g.
#Set pl=T to plot results.  If pl=T, units attribute of y applies.  
#Default is "Day".
#xlab and ... are passed to plot() if pl=T.  Default xlab is label(x)
#if it is defined, otherwise the name of the calling argument for x.
#
#Author: Frank Harrell   8 May 91

groupkm <- function(x, Srv, m=50, g, 
                    cuts, u, pl=FALSE, loglog=FALSE, conf.int=.95, xlab, ylab,
                    lty=1, add=FALSE,
                    cex.subtitle=.7, ...)
{
  if(missing(u))stop("u (time point) must be given")
  if(missing(xlab)) xlab <- label(x)
  if(xlab=="") xlab <- as.character(sys.call())[2]
  s <- !(is.na(x)|is.na(Srv[,1])|is.na(Srv[,2]))
  x <- x[s]; Srv <- Srv[s,]
  x[abs(x) < 1e-10] <- 0 #cut2 doesn't work with tiny x
  e <- Srv[,2]
  if(nrow(Srv) != length(x)) stop("lengths of x and Srv must match")
  unit <- units(Srv)
  if(is.null(unit) || unit=="") unit <- "Day"
  if(!missing(cuts)) q <- cut2(x, cuts)
  else
    if(!missing(g)) q <- cut2(x, g=g)
  else
    q <- cut2(x, m=m)
  if(any(table(q) < 2)) warning('one interval had < 2 observations')

  q <- unclass(q)  #need integer
  g <- length(levels(q))

  km      <- double(g)
  pred    <- km
  std.err <- km
  events  <- integer(g)
  numobs  <- events

#f <- survfit.km(q, Srv, conf.int=conf.int, conf.type="log-log")
#if(is.null(f$strata)) {nstrat <- 1; stemp <- rep(1, length(f$time))}
#else { nstrat <- length(f$strata); stemp <- rep(1:nstrat,f$strata)}
#This is more efficient but doesn't handle empty strata

  for(i in 1:g)
    {
      s <- q==i
      nobs <- sum(s); ne <- sum(e[s])
      if(nobs < 2)
        {
          numobs[i] <- 0
          events[i] <- 0
          pred[i] <- if(nobs==1) mean(x[s], na.rm=TRUE) else NA
          km[i] <- NA
          std.err[i] <- NA
        }
      else
        {
          pred[i] <- mean(x[s], na.rm=TRUE)
          dummystrat <- as.factor(rep(1, nobs))
          f <- survfitKM(dummystrat,Srv[s,]) 
          ##doesn't need conf.int since only need s.e.
          tt <- c(0, f$time)
          ss <- c(1, f$surv)
          se <- c(0, f$std.err)
          tm <- max((1:length(tt))[tt <= u+1e-6])
          km[i] <- ss[tm]
          std.err[i] <- se[tm]
          numobs[i]  <- nobs
          events[i]  <- ne
          n <- length(tt)
          if(u > tt[n]+1e-6 & ss[n]>0)
            {
              km[i] <- NA
              std.err[i] <- NA
            }
        }
    }
  z <- cbind(x=pred, n=numobs, events=events, KM=km, 
             std.err=std.err)

  ciupper <- function(surv, d) ifelse(surv==0, 0, pmin(1, surv*exp(d)))
  cilower <- function(surv, d) ifelse(surv==0, 0, surv*exp(-d))

  if(pl)
    {
      y <- km
      if(conf.int)
        {
          zcrit <- qnorm((conf.int+1)/2)
          low <- cilower(km, zcrit*std.err)
          hi  <- ciupper(km, zcrit*std.err)
        }
      if(missing(ylab))
        ylab <- paste("Kaplan-Meier ",format(u),"-",unit," Survival",sep="")
      if(loglog)
        {
          y <- logb(-logb(y))
          if(conf.int)
            {
              low <- logb(-logb(low))
              hi <- logb(-logb(hi))
            }
          if(missing(ylab))
            ylab <- paste("log(-log Kaplan-Meier ",format(u),unit,
                          " Survival",sep="")
        }
      if(!add)plot(pred, y, xlab=xlab, ylab=ylab, type="n", ...)
      lines(pred, y, lty=lty)
      if(conf.int) errbar(pred, y, hi, low, add=TRUE, ...)
      if(!is.logical(cex.subtitle))
        {
          nn <- sum(numobs,na.rm=TRUE)
          mm <- round(nn/g)
          title(sub=paste("n=",nn," d=",sum(events,na.rm=TRUE),
                  ", avg. ",mm," patients per group",sep=""),
                adj=0,cex=cex.subtitle)
        }
    }
  z
}
