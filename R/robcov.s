robcov <- function(fit, cluster, method=c('huber','efron'))
{
  method <- match.arg(method)

  # Claude Sonnet 5 2026-11-17          5 lines
  # residuals.orm()'s type='score' already correctly uses object$mscore
  # directly (already k+p+qm wide for a yde fit, not just k+p -- see
  # ormll.f90/orm.fit.r) when it is present; it falls back to a
  # different, X-only (p-wide) computation when mscore wasn't
  # requested at fit time, which crashes for a yde fit ("subscript out
  # of bounds") rather than explaining why. Require mscore=TRUE here
  # instead of letting that confusing crash happen -- orm(..., mscore
  # =TRUE) already works today (mscore is forwarded through orm.fit()
  # unchanged, nothing yde-specific needed there).
  if(length(fit$qyde) && fit$qyde > 0 && ! length(fit$mscore))
    stop('robcov() for a yde-fitted orm() object requires the fit to ',
         'have been run with mscore=TRUE')

  var   <- vcov(fit, intercepts='all')
  vname <- dimnames(var)[[1]]
  
  if(inherits(fit, "ols"))
    var <- fit$df.residual * var / sum(fit$residuals ^ 2)  #back to X'X
  else
    if(method=='efron') stop('method="efron" only works for ols fits')

  X <- as.matrix(residuals(fit, type=if(method=='huber')"score" else "hscore"))
  n <- nrow(X)
  if(missing(cluster)) {
    clusterInfo <- NULL
    cluster     <- 1 : n
  }
  else {
    if(any(is.na(cluster))) stop("cluster contains NAs")
    clusterInfo <- list(name=deparse(substitute(cluster)))
  }
  if(length(cluster) != n)
    stop('length of cluster (', length(cluster), ') ',
         'does not match number of observations used in fit (', n, ')')
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
  if(length(clusterInfo)) clusterInfo$n <- length(clus.size)
  clus.start <- c(1, 1 + cumsum(clus.size))
  nc <- length(levels(cluster))
  clus.start <- clus.start[- (nc + 1)]
  storage.mode(clus.start) <- "integer"

  W <- matrix(.Fortran(F_robcovf, n, p, nc, clus.start, clus.size, X, 
                       double(p), w=double(p * p))$w, nrow=p)

##The following has a small bug but comes close to reproducing what robcovf does
##W <- tapply(X,list(cluster[row(X)],col(X)),sum)
##W <- t(W) %*% W
##The following logic will also do it, also at great cost in CPU time
##for(j in levels(cluster))		{
##   s <- cluster==j
##   if(sum(s)==1) sx <- X[s,,drop=F]
##   else {sx <- apply(X[s,,drop=F], 2, sum); dim(sx) <- c(1,p)}
##
##   sc <- sc + t(sx) %*% sx
##
##					}

  adjvar <- var %*% W %*% var
              
  ##var.new <- diag(adjvar)
  ##deff <- var.new/var.orig; names(deff) <- vname
  ##eff.n <- n/exp(mean(log(deff)))

##if(pr)			{
##   v <- cbind(var.orig, var.new, deff)
##   dimnames(v) <- list(vname, c("Original Variance","Adjusted Variance",
##			"Design Effect"))
##   .Options$digits <- 4
##   cat("\n\nEffect of Adjustment for Cluster Sampling on Variances of Parameter #Estimates\n\n")
##   print(v)
##   cat("\nEffective sample size:",format(round(eff.n,1)),"\n\n")
##   nn <- n^2/sum(clus.size^2)
##   cat("\nN^2/[sum of Ni^2]    :",format(round(nn,1)),"\n\n")
##			}

  fit$orig.var <- var
  fit$var <- adjvar
  fit$clusterInfo <- clusterInfo
  ##fit$design.effects <- deff
  ##fit$effective.n <- eff.n
  
  fit
  
}
