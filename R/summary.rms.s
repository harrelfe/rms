
# Value adjusted to is irrelevant when the factor does not interact with
# other factors.  Form of factors is as follows: factor1=value1,factor2=val2:
# Values:
#	NA	: test factor, use all default settings
#	w	: adjust this factor to w when estimating effects of others
#	c(lo,hi): use range for effect (lo,hi), adjust to default value
#	c(lo,w,hi): use range (lo,hi), adjust to w.  Any of 3 can be NA.
# For categories and strata values can be character
# values that are original values before translating to is.category -
# only enough letters are needed to uniquely identify the category 
# This applies to category and strata vars.  Default adjusted to is
# from second element of limits vector.
# For category factors, all comparisons to reference category are made.
# Reference category is assumed to be adjusted to value.
# est.all is T to estimate effects for all factors, not just those listed
# in ...

summary.rms <- function(object, ..., est.all=TRUE, antilog, conf.int=.95,
                           abbrev=FALSE, vnames=c("names","labels"))
{	
  obj.name <- as.character(sys.call())[2]
  at <- object$Design
  labels <- at$label

  vnames <- match.arg(vnames)
  
  assume <- at$assume.code
  if(is.null(assume)) stop("fit does not have design information")
  if(any(assume==10))
    warning("summary.rms does not currently work with matrix factors in model")
  name  <- at$name
  parms <- at$parms

  scale <- object$scale.pred
  if(missing(antilog)) antilog <- length(scale)==2
  if(antilog & length(scale) < 2) scale <- c("","Antilog")

  factors <- rmsArgs(substitute(list(...)))
  nf <- length(factors)

  if(est.all) which <- (1:length(assume))[assume!=9]
  if(nf>0)
    {
      jw <- charmatch(names(factors),name,0)
      if(any(jw==0))stop(paste("factor name(s) not in the design:",
               paste(names(factors)[jw==0],collapse=" ")))
      if(!est.all) which <- jw
      if(any(assume[which]==9))
        stop("cannot estimate effects for interaction terms alone")
    }

  Limval <- Getlim(at, allow.null=TRUE, need.all=FALSE)
  values <- Limval$values
  ## The next statement (9Jun98) makes limits[1:3,] keep all levels of
  ## factors.  Problem is that [.data.frame does not pass drop to []
  ## when first subscripts are specified
  oldopt <- options(drop.factor.levels=FALSE)
  on.exit(options(oldopt))

  lims <- Limval$limits[1:3,,drop=FALSE]
  
  ##Find underlying categorical variables
  ucat <- rep(FALSE, length(assume))
  for(i in (1:length(assume))[assume!=5 & assume<8])
    ucat[i] <- name[i] %in% names(values) &&
  length(V <- values[[name[i]]]) && is.character(V)
  
  stats <- NULL
  lab <- NULL
  lc <- length(object$coef)
  ##Number of non-slopes:
  nrp <- num.intercepts(object)
  nrp1 <- nrp+1
  ## Exclude non slopes
  beta <- object$coef[nrp1:lc]
  var <- vcov(object, regcoef.only=TRUE)[nrp1:lc,nrp1:lc]

  zcrit <- qnorm((1+conf.int)/2)
  cll <- paste(signif(conf.int,3))

  jf <- 0
  if(nf>0) for(i in jw)
    {
      jf <- jf+1
      z <- value.chk(at, i, factors[[jf]], 0, Limval)
      lz <- length(z)
      if(lz==1 && !is.na(z)) lims[2,i] <-  z
      if(lz==2)
        {
          if(!is.na(z[1])) lims[1,i] <- z[1]
          if(!is.na(z[2])) lims[3,i] <- z[2]
        }
      else
        if(lz==3) lims[!is.na(z),i] <- z[!is.na(z)]
      if(lz<1 | lz>3) stop("must specify 1,2, or 3 values for a factor")
    }
  adj <- lims[2,,drop=FALSE]
  isna <- sapply(adj, is.na)


  if(any(isna))
    stop(paste("adjustment values not defined here or with datadist for",
               paste(name[assume!=9][isna],collapse=" ")))
  k <- which[assume[which]!=8 & assume[which]!=5 & assume[which]!=10 & 
             !ucat[which]]
  m <- length(k)
  if(m)
    {
      isna <- is.na(lims[1,name[k],drop=FALSE]+lims[3,name[k],drop=FALSE])
      ##note char. excluded from k
      if(any(isna)) stop(paste("ranges not defined here or with datadist for",
                               paste(name[k[isna]], collapse=" ")))
    }
  
  xadj <- oldUnclass(rms.levels(adj, at))
  m <- length(k)
  if(m)
    {
      adj <- xadj
      M <- 2*m
      odd <- seq(1,M,by=2)
      even<- seq(2,M,by=2)
      ##Extend data frame
      for(i in 1:length(adj)) adj[[i]] <- rep(adj[[i]], M)
   
      i <- 0
      for(l in k)
        {
          i <- i+1
          adj[[name[l]]][(2*i-1):(2*i)] <- lims[c(1,3),name[l]]
        }
      xx <- predictrms(object, newdata=adj, type="x", incl.non.slopes=FALSE)
      xd <- matrix(xx[even,]-xx[odd,],nrow=m)
      xb <- xd %*% beta
      se <- drop((((xd %*% var) * xd) %*% rep(1,ncol(xd)))^.5)
      low <- xb - zcrit*se
      up <- xb + zcrit*se
      lm <- as.matrix(lims[,name[k],drop=FALSE])
      stats <- cbind(lm[1,],lm[3,],lm[3,]-lm[1,],xb,se,low,up,1)
      lab <- if(vnames=='names') name[k] else labels[k]
      if(antilog)
        {
          stats <- rbind(stats,
                         cbind(stats[,1:3,drop=FALSE],
                               exp(xb),NA,exp(low),exp(up), 2))
          lab <- c(lab,rep(paste("",scale[2]),m))
          w <- integer(M)
          w[odd] <- 1:m
          w[even]<- m+(1:m)
          stats <- stats[w,]
          lab <- lab[w]
        }
    }
  
  for(j in 1:length(xadj)) xadj[[j]] <- rep(xadj[[j]], 2)

  for(i in which[assume[which]==5 | ucat[which]])
    {
      ## All comparisons with reference category
  
      parmi <- if(ucat[i]) values[[name[i]]] else parms[[name[i]]]
      parmi.a <- if(abbrev) abbreviate(parmi) else parmi
      iref <- as.character(xadj[[name[i]]][1])
      ki <- match(iref, parmi)
      for(j in parmi)
        {
          if(j!=iref)
            {
              kj <- match(j, parmi)
              adj <- xadj
              adj[[name[i]]] <- c(iref,j)
              adj <- as.data.frame(adj)
              xx <- predictrms(object,newdata=adj,
                               type="x",incl.non.slopes=FALSE)
              xd <- matrix(xx[2,]-xx[1,],nrow=1)
              xb <- (xd %*% beta)
              se <- sqrt((xd %*% var) %*% t(xd))
              low <- xb - zcrit*se
              up <- xb + zcrit*se
              stats <- rbind(stats,cbind(ki,kj,NA,
                                         xb,se,low,up,1))
              lab <-c(lab,
                      paste(if(vnames=='names') name[i] else labels[i],
                            " - ",parmi.a[kj],":",
                            parmi.a[ki],sep=""))
              if(antilog)
                {
                  stats <- rbind(stats,cbind(ki,kj,NA,
                                             exp(xb),NA,exp(low),exp(up),2))
                  lab <- c(lab, paste("",scale[2]))}
            }
        }
    }

  dimnames(stats) <-
    list(lab, c("Low","High",
                "Diff.","Effect","S.E.",
                paste("Lower",cll),paste("Upper",cll),"Type"))
  
  attr(stats,"heading") <-
    paste("             Effects              Response : ",
          as.character(formula(object))[2], sep='')
  attr(stats,"class") <- c("summary.rms","matrix")
  attr(stats,"scale") <- scale
  attr(stats,"obj.name") <- obj.name
  interact <- at$interactions
  adjust <- ""
  if(length(interact))
    {
      interact <- sort(unique(interact[interact>0]))
      nam <- name[which[match(which, interact, 0)>0]]
      if(length(nam)) for(nm in nam) 
        adjust <- paste(adjust, nm,"=",
                        if(is.factor(xadj[[nm]]))
                        as.character(xadj[[nm]])[1] else
                        format(xadj[[nm]][1])," ",sep="")
    }
  attr(stats,"adjust") <- adjust
  
  stats
}



print.summary.rms <- function(x, ...)
{
  cstats <- dimnames(x)[[1]]
  for(i in 1:3) cstats <- cbind(cstats, format(signif(x[,i],5)))
  for(i in 4:7) cstats <- cbind(cstats, format(round(x[,i],2)))
  dimnames(cstats) <- list(rep("",nrow(cstats)), 
                           c("Factor", dimnames(x)[[2]][1:7]))
  cat(attr(x,"heading"),"\n\n")
  print(cstats,quote=FALSE)
  if((A <- attr(x,"adjust"))!="") cat("\nAdjusted to:", A,"\n\n")
  
  invisible()
}


latex.summary.rms <-
  function(object, 
           title=if(under.unix)
           paste('summary',attr(object,'obj.name'),sep='.') else
           paste("sum",substring(first.word(attr(object,"obj.name")),
                                 1,5),sep=""),
           table.env=TRUE, ...)
{ 

  title <- title   # because of lazy evaluation
  caption <- latexTranslate(attr(object, "heading"))
  scale <- attr(object,"scale")
  object <- object[,-8,drop=FALSE]
  rowl <- latexTranslate(dimnames(object)[[1]])
  rowl <- ifelse(substring(rowl,1,1)==" ",
                 paste("~~{\\it ",
                       substring(rowl,2),"}",sep=""),
                 rowl) # preserve leading blank
  rowl <- sedit(rowl, "-", "---")
  cstats <- matrix("", nrow=nrow(object), ncol=ncol(object), 
                   dimnames=dimnames(object))
  for(i in 1:3) cstats[,i] <- format(signif(object[,i],5))
  for(i in 4:7) cstats[,i] <- format(round(object[,i],2))
  cstats[is.na(object)] <- ""
  caption <- sedit(caption, "    Response","~~~~~~Response")
  cstats <- as.data.frame(cstats)
  attr(cstats,"row.names") <- rowl
  names(cstats)[3] <- "$\\Delta$"
  latex(cstats, caption=if(table.env) caption else NULL,
        title=title, rowlabel="",
        col.just=rep("r",7), table.env=table.env, ...)
}


plot.summary.rms <-
  function(x, at, log=FALSE, 
           q=c(0.7, 0.8, 0.9, 0.95, 0.99), xlim, nbar, cex=1, nint=10, cex.c=.5,
           cex.t=1, clip=c(-1e30,1e30), main, ...)
{
  scale  <- attr(x, "scale")
  adjust <- attr(x, "adjust")

  Type   <- x[,"Type"]
  x  <- x[Type==1,,drop=FALSE]
  lab    <- dimnames(x)[[1]]
  effect <- x[,"Effect"]
  se     <- x[,"S.E."]
  if(!log && any(Type==2))
    {
      fun <- exp
      tlab <- scale[2]
    }
  else
    {
      fun <- function(x) x
      if(log)
        {
          if(length(scale)==2) tlab <- scale[2]
          else tlab <- paste("exp(",scale[1],")",sep="")
        }
      else tlab <- scale[1]
    }
  if(!length(scale)) tlab <- ''  ## mainly for Glm fits
  if(!missing(main)) tlab <- main
  augment <- if(log | any(Type==2)) c(.1, .5, .75, 1) else 0
  n     <- length(effect)
  out   <- qnorm((max(q)+1)/2)
  if(missing(xlim) && !missing(at)) xlim <- range(if(log)logb(at) else at) else
  if(missing(xlim))
    {
      xlim <- fun(range(c(effect-out*se,effect+out*se)))
      xlim[1] <- max(xlim[1],clip[1])
      xlim[2] <- min(xlim[2],clip[2])
    }
  else
    augment <- c(augment, if(log)exp(xlim) else xlim)
  
  fmt <- function(k)
    {
      m <- length(k)
      f <- character(m)
      for(i in 1:m) f[i] <- format(k[i])
      f
    }
  lb <- ifelse(is.na(x[,'Diff.']), lab,
               paste(lab,' - ',
                     fmt(x[,'High']),':',fmt(x[,'Low']),sep=''))
  plot.new(); par(new=TRUE)
  mxlb <- .1+max(strwidth(lb,units='inches',cex=cex))
  tmai <- par('mai')
  on.exit(par(mai=tmai))
  par(mai=c(tmai[1],mxlb,1.5*tmai[3],tmai[4]))
  
  outer.widths <- fun(effect+out*se)-fun(effect-out*se)
  if(missing(nbar)) nbar <- n
  npage <- ceiling(n/nbar)
  is <- 1
  for(p in 1:npage)
    {
      ie <- min(is+nbar-1, n)
      plot(1:nbar, rep(0,nbar), xlim=xlim, ylim=c(1,nbar),
           type="n", axes=FALSE, 
           xlab="", ylab="")
      if(cex.t>0) title(tlab, cex=cex.t)
      lines(fun(c(0,0)),c(nbar-(ie-is), nbar),lty=2)
      if(log)
        {
          pxlim <- pretty(exp(xlim), n=nint)
          pxlim <- sort(unique(c(pxlim, augment)))
          ## For wome weird reason, sometimes duplicates (at xlim[2])
          ## still remain
          pxlim <- pxlim[pxlim>=exp(xlim[1])]
          if(!missing(at)) pxlim <- at
          axis(3, logb(pxlim), labels=format(pxlim))
        }
      else 
        {
          pxlim <- pretty(xlim, n=nint)
          pxlim <- sort(unique(c(pxlim, augment)))
          pxlim <- pxlim[pxlim>=xlim[1]]
          if(!missing(at)) pxlim <- at
          axis(3, pxlim)
        }
      imax <- (is:ie)[outer.widths[is:ie]==max(outer.widths[is:ie])][1]
      for(i in is:ie)
        {
          confbar(nbar-(i-is+1)+1, effect[i], se[i], q=q, type="h", 
                  fun=fun, cex=cex.c, labels=i==imax, clip=clip, ...)
          mtext(lb[i], 2, 0, at=nbar-(i-is+1)+1, cex=cex,
                adj=1, las=1)
        }
      if(adjust!="") 
        {
          adjto <- paste("Adjusted to:",adjust,sep="")
          xx <- par('usr')[2]
          if(nbar>ie) text(xx, nbar-(ie-is+1), adjto, adj=1, cex=cex)
          else title(sub=adjto, adj=1, cex=cex)
        }
      is <- ie+1
    }
  invisible()
}
