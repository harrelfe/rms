### The idea is to get score residuals from the score vector ###

# a function to get score residiuals from the fitted model
# (the code are taking directly from `orm.fit`)
func_score <- function(fit){
  
  # coefficients
  coef <- fit$coefficients
  
  # convert y to ordered category
  y <- match(fit$y, fit$yunique)
  # x
  x <- fit$x
  
  # some useful numbers
  kint <- length(fit$yunique) - 1L
  nx <- dim(x)[2L]
  n <- length(y)
  p <- as.integer(kint + nx)
  
  # store results - matrix
  u <- matrix(0, nrow = n, ncol=p)
  
  # functions
  f   <- eval(fit$famfunctions[1])
  fp  <- eval(fit$famfunctions[3])
  
  
  xb <- fit$x %*% coef[-(1L : kint)]
  ints <- c(1e100, coef[1:kint], -1e100)
  xby <- xb + ints[y]; xby1 <- xb + ints[y + 1L]
  fa <- f(xby)
  fb <- f(xby1)
  P <- fa - fb
  fpa  <- fp(xby,   fa)
  fpb  <- fp(xby1,  fb)
  
  # score for alpha
  for(m in 1:kint){
    for(j in 1:n){
      u[j, m] <-  (fpa[j]*(y[j]-1==m) - fpb[j]*(y[j]==m)) / P[j]
    }
  }
  
  # score for beta
  for(m in (kint+1):p){
    for(j in 1:n){
      u[j, m] <- (fpa[j] - fpb[j]) * x[j,m-kint] / P[j]
    }
  }
  
  return(u)
  
}


func_robcov <- function(fit, cluster){
  var   <- vcov(fit, intercepts='all')
  vname <- dimnames(var)[[1]]
  
  # X <- as.matrix(residuals(fit, type="score"))
  X <- func_score(fit) # get score residuals
  
  n <- nrow(X)
  
  cluster <- as.factor(cluster)
  
  p <- ncol(var)
  j <- is.na(X %*% rep(1, ncol(X)))
  if(any(j)) {
    X       <- X[! j,, drop=FALSE]
    cluster <- cluster[! j, drop=TRUE]
    n       <- length(cluster)
  }
  
  j <- order(cluster)
  X <- X[j, , drop=FALSE]
  clus.size  <- table(cluster)
  # if(length(clusterInfo)) clusterInfo$n <- length(clus.size)
  clus.start <- c(1, 1 + cumsum(clus.size))
  nc <- length(levels(cluster))
  clus.start <- clus.start[- (nc + 1)]
  storage.mode(clus.start) <- "integer"
  
  # dyn.load("robcovf.so")
  W <- matrix(.Fortran("robcovf", n, p, nc, clus.start, clus.size, X, 
                       double(p), w=double(p * p))$w, nrow=p)
  
  ##The following has a small bug but comes close to reproducing what robcovf does
  # W <- tapply(X,list(cluster[row(X)],col(X)),sum)
  # W <- t(W) %*% W
  
  #The following logic will also do it, also at great cost in CPU time
  # W <- matrix(0, p, p)
  # for(j in levels(cluster)){
  #   s <- cluster==j
  #   if(sum(s)==1){
  #     sx <- X[s,,drop=F]
  #   }else {sx <- apply(X[s,,drop=F], 2, sum); dim(sx) <- c(1,p)}
  #   W <- W + t(sx) %*% sx
  # }
  
  adjvar <- var %*% W %*% var
  
  return(adjvar)
}



### test ###
if(FALSE) {
# generate longitudinal data
library(mvtnorm)
library(rms)
data_gen <- function(n=100, m=10, b=beta, a0=0.7){
  # correlation matrix - exchangeable structure
  G <- matrix(rep(a0, m*m), nrow=m)
  diag(G) <- 1
  
  stdevs <- rep(1,m)
  e <- rmvnorm(n, mean = rep(0,m), 
               sigma = G * matrix(outer(stdevs, stdevs), nrow=m, byrow=TRUE))

  # x1: gender (0: female, 1: male)
  x1 <- rep(c(rep(0,round(n/2)), rep(1,n-round(n/2))), m)
  
  # x2: time
  x2 <- rep(c(1:m), each=n)
  
  y <- b[1] + b[2] * x1 + b[3] * x2 + e
  
  dat <- data.frame(y=c(y), 
                    x1=c(x1), 
                    x2=c(x2), 
                    id=rep(1:n, m))
  
  return(dat)
}

# data
dat <- data_gen(n=50, m=10, b=beta, a0=0.7)
# model
mod_orm <- orm(y ~ x1 + x2,
               data = dat,
               x=T, y=T)
# robcov
rob.cov <- robcov(fit = mod_orm,
                  cluster = dat$id)
rob.cov
}
