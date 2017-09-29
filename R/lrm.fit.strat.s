lrm.fit.strat <- function(x, y, strata, offset=0, initial,
						  maxit=25, eps=.025, tol=1E-7, trace=FALSE,
						  penalty.matrix=NULL, strata.penalty=0,
                          weights=NULL, normwt) {
  cal <- match.call()
  opts <- double(12)
  len.penmat <- length(penalty.matrix)
  lev    <- levels(strata)
  nstrat <- length(lev)
  strata <- unclass(strata)
  n <- length(y)
  
  if(!length(weights)) {
    normwt <- FALSE
    weights <- rep(1,n)
  }
  if(length(weights) != n) stop('weights and y must be same length')
  storage.mode(weights) <- 'double'
  opts[12] <- normwt
  ## weights not implemented for stratified models yet

  initial.there <- !missing(initial)
  if(missing(x) || length(x)==0) {
    nx <- 0
    xname <- NULL
    x <- 0
  }
  else {
    if(!is.matrix(x)) x <- as.matrix(x)
    storage.mode(x) <- "double"
    dx <- dim(x)
    nx <- dx[2]
    if(dx[1]!=n)stop("x and y must have same length")
    xname <- dimnames(x)[[2]]
    if(length(xname)==0) xname <- paste("x[",1:nx,"]",sep="")
  }
  
  nxin <- nx
  
  if(!is.factor(y)) y <- as.factor(y)
  y <- unclass(y)   # in case is.factor
  ylevels <- levels(y)
  
  if(n < 3)stop("must have >=3 non-missing observations")
  kint <- as.integer(length(ylevels)-1)
  if(kint != 1) stop('only works for binary y')
  ftable <- integer(5001*(kint+1))
  levels(y) <- ylevels
  numy <- table(y)
  y <- as.integer(y-1)
  nvi <- as.integer(nxin+kint+nstrat-1)
  if(missing(initial)) {
    ncum <- rev(cumsum(rev(numy)))[2:(kint+1)]
    pp <- ncum/n
    initial <-logb(pp/(1-pp))
    initial <- initial - mean(offset)
  }
  if(length(initial) < nvi) initial <- c(initial,rep(0,nvi-length(initial)))
  storage.mode(initial) <- "double"
  loglik <- -2 * sum(numy * logb(numy/n))
  
  if(nxin > 0) {
    if(len.penmat==0) penalty.matrix <- matrix(0,nrow=nx,ncol=nx)
    if(nrow(penalty.matrix)!=nx || ncol(penalty.matrix)!=nx) 
      stop(paste("penalty.matrix does not have",nx,"rows and columns"))
    penmat <- rbind(
                matrix(0,ncol=kint+nx,nrow=kint),
                cbind(matrix(0,ncol=kint,nrow=nx),penalty.matrix))
  }
  else
    penmat <- matrix(0, ncol=kint, nrow=kint)
  storage.mode(penmat) <- 'double'
  
  ofpres <- !all(offset == 0)
  storage.mode(offset) <- 'double'
  
  if(nxin==0 & !ofpres)
    {
      loglik <- rep(loglik,2)
      z <- list(coef=initial,u=rep(0,kint),opts=c(rep(0,7),.5,0,0,0))
    }

  if(ofpres)
    {
      ##Fit model with only intercept(s) and offset
      z <- 
        .Fortran("lrmfit",coef=initial,as.integer(0),0,x,y,offset,
                 u=double(kint),
                 double(kint*(kint+1)/2),loglik=double(1),n,as.integer(0),
                 numy,kint,
                 v=double(kint*kint),double(kint),double(kint),
                 double(kint),pivot=integer(kint),opts=opts,ftable,
                 penmat,weights, PACKAGE="rms")
      
      loglik <- c(loglik,z$loglik)
      if(z$opts[6] | z$opts[7]<kint) return(list(fail=TRUE,class="lrm"))
      initial <- z$coef
    }

		
  ## Newton-Raphson iterations with patterned matrix inversion for speed
  theta  <- initial
  iter   <- 0
  oldobj <- 1e10
  x      <- cbind(1,x)
  nns    <- nx+1  ## no. of non-strata parameters
  
  while(iter <= maxit)
    {
      iter <- iter + 1
      beta <- as.matrix(theta[1:nns])
      tau  <- c(0,theta[-(1:nns)])
      logit <- drop(x %*% beta + tau[strata] + offset)
      pred <- 1/(1+exp(-logit))
      obj  <- -2*sum(y*logb(pred) + (1-y)*logb(1-pred)) + 
        t(beta) %*% penmat %*% beta + 
          strata.penalty*sum((tau-mean(tau))^2)
      if(trace)cat('-2logL=',format(obj),'')
      d <- y - pred
      u <-  c(drop(matrix(d,nrow=1) %*% x)-drop(penmat %*% beta), 
              tapply(d,strata,sum)[-1] - strata.penalty*(tau[-1]-mean(tau)))
      u <- as.matrix(u)
      if(trace)cat(format(u),'\n')
      ## Created patterned information matrix  A B / B' C
      ## Inverse is AA BB / BB' CC
      pq <- pred*(1-pred)
      A <- crossprod(pq * x, x) + penmat
      B <- t(rowsum(pq * x, strata))[,-1,drop=FALSE]  
      ## above won't work if a stratum not represented
      dd <- tapply(pq, strata, sum)[-1]
      v <- 1/(dd + strata.penalty)
      vm <- as.matrix(v)
      k <- (strata.penalty/nstrat)/(1 - (strata.penalty/nstrat)*sum(v))
      ## BCi <- t(v * t(B)) + k * (B %*% vm) %*% t(vm)
      ## C <- diag(dd + strata.penalty) - strata.penalty/nstrat
      BCi <- B*rep(v,rep(nns,nstrat-1)) + k * (B %*% vm) %*% t(vm)
      AA <- solvet(A - BCi %*% t(B), tol=tol)
      BB <- -AA %*% BCi
      u1 <- u[1:nns,,drop=FALSE]
      u2 <- u[(nns+1):nvi,,drop=FALSE]
      theta <- theta + c(AA %*% u1 + BB %*% u2,
                         t(BB) %*% u1 + vm * u2 + k * vm %*% (t(vm) %*% u2) -
                         t(BB) %*% (BCi %*% u2))
      ##CC <- diag(drop(vm))+k*vm%*%t(vm)-t(BB) %*% BCi
      ##					 t(BB) %*% u1 + u2/dd)             FAILS
      ##theta <- theta + c(AA %*% u1 + BB %*% u2,
      ##					 t(BB) %*% u1 + (diag(1/dd) - t(BB) %*% BCi) %*% u2)
      ## theta <- theta + c(solve(A) %*% u1, u2/dd)           SLOW
      ##  theta <- theta + c(diag(1/diag(A)) %*% u1, u2/dd)   FAILS
      if(abs(obj - oldobj) < eps) break
      oldobj <- obj
    }
  if(iter > maxit) return(list(fail=TRUE, class='lrm'))

  
  xname <- c(xname, lev[-1])

  if(kint==1) name <- "Intercept"
  else 
    name <- paste("y>=",ylevels[2:(kint+1)],sep="")
  name <- c(name, xname)
  theta <- drop(theta)
  names(theta) <- name
  
  loglik <- c(loglik, obj)
  
  dimnames(AA) <- list(name[1:nns],name[1:nns])
  dimnames(BB) <- dimnames(BCi) <- list(name[1:nns],name[(nns+1):nvi])
  names(BCi)   <- NULL
  
  
  llnull <- loglik[length(loglik)-1]
  model.lr <- llnull-loglik[length(loglik)]
  model.df <- nvi - kint
  if(initial.there) model.p <- NA
  else
    {
      if(model.df>0) model.p <- 1-pchisq(model.lr,model.df)
      else
        model.p <- 1
    }
  r2 <- 1-exp(-model.lr/n)
  r2.max <- 1-exp(-llnull/n)
  r2 <- r2/r2.max
  Brier <- mean((pred - (y>0))^2)
  
  stats <- c(n,max(abs(u)),model.lr,model.df,model.p,
             ## z$opts[8],z$opts[9],z$opts[10], z$opts[11], 
             r2, Brier)
  nam <- c("Obs","Max Deriv",	"Model L.R.","d.f.","P",
           ##"C","Dxy","Gamma","Tau-a",
           "R2","Brier")
  names(stats) <- nam

  vcov <- function(fit, which=c('strata.var','var','strata.var.diag'))
    {
      which <- match.arg(which)
      strata.penalty <- fit$strata.penalty
      v <- 1 / (fit$strata.unpen.diag.info + strata.penalty)
      nstrat <- fit$nstrat
      k <- (strata.penalty/nstrat)/(1 - (strata.penalty/nstrat)*sum(v))
      sname <- fit$strata.levels[-1]
      CC <- diag(v) + k * v %*% t(v) -t(fit$cov.nonstrata.strata) %*% fit$BCi
      switch(which,
             strata.var = structure(CC, dimnames=list(sname,sname)),
             strata.var.diag = structure(diag(CC), names=sname),
             var = structure(rbind(cbind(fit$var,fit$cov.nonstrata.strata),
               cbind(t(fit$cov.nonstrata.strata),CC)),
               dimnames=list(nn <- names(fit$coef),nn)))
    }

  retlist <- list(call=cal,freq=numy,
                  stats=stats,fail=FALSE,coefficients=theta[1:nns],
                  non.slopes=1,est=1:(nvi-kint),
                  var=AA,u=u,
                  deviance=loglik,
                  linear.predictors=logit,
                  penalty.matrix=if(nxin>0 && any(penalty.matrix!=0)) 
				  penalty.matrix else NULL,
                  nstrat=nstrat, strata.levels=lev,
                  strata.coefficients=theta[(nns+1):nvi],
                  strata.penalty=strata.penalty, 
                  strata.unpen.diag.info=dd,
                  cov.nonstrata.strata=BB,
                  BCi=BCi,
                  vcov=vcov,
                  ## info.matrix=rbind(cbind(A,B),cbind(t(B),diag(dd))))
                  info.matrix=A)
  class(retlist) <- c("lrm","lm")
  retlist
}

