sensuc <- function(fit,
				   or.xu=seq(1,6,by=.5), or.u=or.xu, prev.u=.5,
				   constrain.binary.sample=TRUE,
				   or.method=c('x:u y:u','u|x,y'),
				   event=function(y) if(is.matrix(y))y[,ncol(y)] else 1*y)
{
  type <- class(fit)[1]
  if(type %nin% c('lrm','cph')) stop('fit must be from lrm or cph')
  
  or.method <- match.arg(or.method)
  
  X <- fit$x
  Y <- fit$y
  if(length(X)==0 || length(Y)==0) stop('did not specify x=TRUE, y=TRUE to fit')
  x <- X[,1]
  unq <- sort(unique(x))
  if(length(unq) != 2 || unq[1] != 0 || unq[2] != 1)
	stop('x is not binary')
  
  
  event <- event(Y)
  unq <- sort(unique(event))
  if(length(unq) != 2 || unq[1] != 0 || unq[2] != 1)
	stop('Y or event is not binary')
  
  ##Function to generate Bernoullis with exact proportion p except for roundoff
  bern <- function(n, p, constrain)
    {
      if(constrain)
        {
          sort.random <- function(x)
            {
              un <- runif(length(x))
              x[order(un)]
            }
          ones <- round(n*p)
          zeros <- n - ones
          sort.random(c(rep(0,zeros),rep(1,ones)))
        }
      else
        sample(0:1, n, replace=TRUE, c(1-p,p))
    }

  a00 <- mean(!event & !x)
  a10 <- mean(event & !x)
  a01 <- mean(!event & x)
  a11 <- mean(event & x)
  p.event <- mean(event)
  
  b1  <- p.event
  b0  <- 1 - b1
  c1  <- mean(x)
  c0  <- 1 - c1
  
  n <- length(event)
  
  n00 <- sum(!event & !x)
  n10 <- sum(event & !x)
  n01 <- sum(!event & x)
  n11 <- sum(event & x)
  m1  <- prev.u * n
  m0  <- n - m1



  m <- length(or.xu) * length(or.u)

  OR.xu <- OR.u <- effect.x <- OOR.xu <- effect.u <- effect.u.adj <- Z <- 
	double(m)

  Prev.u <- matrix(NA,nrow=m,ncol=4,
				   dimnames=list(NULL,c('event=0 x=0','event=1 x=0',
					 'event=0 x=1','event=1 x=1')))

  odds <- function(x)
    {
      p <- mean(x)
      p/(1-p)
    }
  
  j <- 0
  cat('Current odds ratio for x:u=')
  for(c.or.xu in or.xu)
    {
      cat(c.or.xu,'')
      for(c.or.u in or.u)
        {
          j <- j + 1
          OR.xu[j] <- c.or.xu
          OR.u[j]  <- c.or.u
          
          if(or.method=='u|x,y')
            {
              beta  <- logb(c.or.u)
              gamma <- logb(c.or.xu)
              f <- function(alpha,beta,gamma,a00,a10,a01,a11,prev.u) 
                a00*plogis(alpha)+
                  a10*plogis(alpha+beta)+
                    a01*plogis(alpha+gamma)+
                      a11*plogis(alpha+beta+gamma) - prev.u
              
              alpha <- uniroot(f, lower=-10, upper=10,
                               beta=beta, gamma=gamma, 
                               a00=a00, a10=a10, a01=a01, a11=a11, 
                               prev.u=prev.u)$root
              p00 <- plogis(alpha)
              p10 <- plogis(alpha+beta)
              p01 <- plogis(alpha+gamma)
              p11 <- plogis(alpha+beta+gamma)
            }
          else
            {
              ## Raking method, thanks to M Conaway
              rake2x2 <- function(prow,pcol,odds) {
                pstart <- matrix(1, nrow=2, ncol=2)
                pstart[1,1] <- odds
                pstart <- pstart/sum(pstart)
                oldp <- pstart
                maxdif <- 1
                while(maxdif > .0001)
                  {
                    ## Adjust row totals
                    obsrow <- oldp[,1]+oldp[,2]
                    adjrow <- prow / obsrow
                    newp <- oldp * cbind(adjrow,adjrow)
                    ## Adjust col totals
                    obscol <- newp[1,]+newp[2,]
                    adjcol <- pcol / obscol
                    newp <- newp * rbind(adjcol,adjcol)
                    maxdif <- max(abs(newp - oldp))
                    oldp <- newp
                  }
                c(newp[1,],newp[2,])
              }
		
              lambda <- c.or.xu
              theta  <- c.or.u
              prow <- c(1-prev.u, prev.u)
              pcol <- c(n00,n01,n10,n11)/n
              a <- matrix(c(
                            1,0,1,0,0,0,0,0,
                            0,1,0,1,0,0,0,0,
                            0,0,0,0,1,0,1,0,
                            0,0,0,0,0,1,0,1,
                            1,1,0,0,0,0,0,0,
                            0,0,1,1,0,0,0,0,
                            0,0,0,0,1,1,0,0,
                            0,0,0,0,0,0,1,1,
                            1,0,0,0,1,0,0,0,
                            0,1,0,0,0,1,0,0,
                            0,0,1,0,0,0,1,0,
                            0,0,0,1,0,0,0,1),
                          nrow=12,byrow=TRUE)
              aindx <- matrix(c(
                                1,3,
                                2,4,
                                5,7,
                                6,8,
                                1,2,
                                3,4,
                                5,6,
                                7,8,
                                1,5,
                                2,6,
                                3,7,
                                4,8),
                              ncol=2, byrow=TRUE)
              pcol1 <- c(pcol[1]+pcol[3], pcol[2]+pcol[4])
              u <- rake2x2(prow, pcol1, lambda)
              
              pcol2 <- c(pcol[1]+pcol[2],pcol[3]+pcol[4])
              w <- rake2x2(prow, pcol2, theta)
              
              newp8 <- p8start <- rep(1/8, 8)
              targvec <- c(u, w, pcol)
              d <- 1
              while(d > .0001)
                {
                  for(i in 1:12)
                    {
                      adjust <- targvec[i] / sum(a[i,] * newp8)
                      newp8[aindx[i,]] <- adjust * newp8[aindx[i,]]
                    }
                  chktarg <- a %*% as.matrix(newp8)
                  d <- max(abs(chktarg - targvec))
                }
              p00 <- newp8[5]/a00
              p01 <- newp8[6]/a01
              p10 <- newp8[7]/a10
              p11 <- newp8[8]/a11
              
##		  prn(c(newp8[5],newp8[5]*n,newp8[5]/(newp8[1]+newp8[5]),
##				newp8[5]*n/n00,newp8[5]/a00))
##		  w_newp8
##		  A_w[1];B_w[2];C_w[3];D_w[4];E_w[5];FF_w[6];G_w[7];H_w[8]
##		  prn((FF+H)*(A+C)/(B+D)/(E+G))
##		  prn((G+H)*(A+B)/(E+FF)/(C+D))
##		  w1_p01*b0+p11*b1
##		  w2_p00*b0+p10*b1
##		  prn((w1/(1-w1))/(w2/(1-w2)))
##		  z1_p10*c0+p11*c1
##		  z2_p00*c0+p10*c1
##		  prn((z1/(1-z1))/(z2/(1-z2)))
            }
	  
          Prev.u[j,] <- c(p00,p10,p01,p11)

          u <- rep(0, n)
          i <- !event & !x
          u[i] <- bern(sum(i), p00, constrain.binary.sample)
          i <- event & !x
          u[i] <- bern(sum(i), p10, constrain.binary.sample)
          i <- !event & x
          u[i] <- bern(sum(i), p01, constrain.binary.sample)
          i <- event & x
          u[i] <- bern(sum(i), p11, constrain.binary.sample)
          
          OOR.xu[j] <-  odds(u[x==1])/odds(u[x==0])
          
          if(type=='cph')
            {
              g <- coxphFit(as.matrix(u),Y,rep(1,n),toler.chol=1e-11,
                            iter.max=15,eps=.0001,method='efron')
              effect.u[j] <- exp(g$coefficients)
              
              g <- coxphFit(cbind(u,X),Y,rep(1,n),toler.chol=1e-11,
                            iter.max=15,eps=.0001,method='efron')
              cof <- g$coefficients
              vr  <- g$var
            }
          else
            {
              effect.u[j] <- odds(event[u==1])/odds(event[u==0])
              g <- lrm.fit(cbind(u,X),event,maxit=20,tol=1E-11)
              ns <- g$non.slopes
              cof <- g$coefficients[-(1:ns)]
              vr  <- g$var[-(1:ns),-(1:ns)]
            }
          z   <- cof/sqrt(diag(vr))
          
          effect.u.adj[j] <- exp(cof[1])
          effect.x[j] <- exp(cof[2])
          Z[j]    <- z[2]
        }
    }
  cat('\n\n')
  
  structure(list(OR.xu=OR.xu,OOR.xu=OOR.xu,OR.u=OR.u,
				 effect.x=effect.x,effect.u=effect.u,effect.u.adj=effect.u.adj,
				 Z=Z,prev.u=prev.u,cond.prev.u=Prev.u,
				 type=type), class='sensuc')
}

plot.sensuc <- function(x, ylim=c((1+trunc(min(x$effect.u)-.01))/
                             ifelse(type=='numbers',2,1),
                             1+trunc(max(x$effect.u)-.01)),
						xlab='Odds Ratio for X:U',
						ylab=if(x$type=='lrm')'Odds Ratio for Y:U' else
						'Hazard Ratio for Y:U',
						digits=2, cex.effect=.75, cex.z=.6*cex.effect,
						delta=diff(par('usr')[3:4])/40, 
						type=c('symbols','numbers','colors'),
						pch=c(15,18,5,0),col=c(2,3,1,4),alpha=.05,
						impressive.effect=function(x)x > 1,...)
{
  type  <- match.arg(type)
  Z     <- abs(x$Z)
  or    <- x$OOR.xu
  eu    <- x$effect.u
  ex    <- x$effect.x
  zcrit <- qnorm(1-alpha/2)

  plot(or, eu, ylim=ylim, xlab=xlab, ylab=ylab, type='n', ...)

  if(type=='numbers')
    {
      text(or, eu, round(ex,digits), cex=cex.effect)
      text(or, eu - delta, round(Z,2), cex=cex.z)
    }
  else
    {
      i <- impressive.effect(ex) & Z >= zcrit
      if(any(i)) if(type=='symbols') points(or[i], eu[i], pch=pch[1])
      else
        text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[1])
      i <- impressive.effect(ex) & Z < zcrit
      if(any(i))
        if(type=='symbols') 
          points(or[i], eu[i], pch=pch[2])
        else
          text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[2])
      i <- !impressive.effect(ex) & Z < zcrit
      if(any(i)) if(type=='symbols') points(or[i], eu[i], pch=pch[3])
      else
        text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[3])
      i <- !impressive.effect(ex) & Z >= zcrit
      if(any(i))
        if(type=='symbols')
          points(or[i], eu[i], pch=pch[4])
        else
          text(or[i], eu[i], round(ex[i],digits), cex=cex.effect, col=col[4])
    }
  title(sub=paste('Prevalence of U:',format(x$prev.u)),adj=0)
  invisible()
}
