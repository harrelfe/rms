##' Produce Design Matrices for Contrasts
##'
##' This is a simpler version of `contrast.rms` that creates design matrices or differences of them and does not require the fit object to be complete (i.e., to have coefficients).  This is used for the `pcontrast` option in [rmsb::blrm()].
##' @title Xcontrast
##' @param fit an `rms` or `rmsb` fit object, not necessarily complete
##' @param a see [rms::contrast.rms()]
##' @param b see [rms::contrast.rms()]
##' @param a2 see [rms::contrast.rms()]
##' @param b2 see [rms::contrast.rms()]
##' @param ycut see [rms::contrast.rms()]
##' @param weights see [rms::contrast.rms()]
##' @param expand see [rms::contrast.rms()]
##' @return numeric matrix
##' @author Frank Harrell
Xcontrast <- function(fit, a, b=NULL, a2=NULL, b2=NULL, ycut=NULL,
                      weights='equal', expand=TRUE) {
    
  partialpo <- inherits(fit, 'blrm') && fit$pppo > 0
  if(partialpo & ! length(ycut))
    stop('must specify ycut for partial prop. odds model')
  cppo      <- fit$cppo
  if(partialpo && ! length(cppo))
    stop('only implemented for constrained partial PO models')
  
  pred <- function(d) {
    ## predict.blrm duplicates rows of design matrix for partial PO models
    ## if ycut has length > 1 and only one observation is being predicted
    if(partialpo) predict(fit, d, type='x', ycut=ycut)
         else
           predict(fit, d, type='x')
    }
  
  da <- do.call('gendata', list(fit, factors=a, expand=expand))
  xa <- pred(da)
  if(length(b)) {
    db <- do.call('gendata', list(fit, factors=b, expand=expand))
    xb <- pred(db)
    }

  ma <- nrow(xa)

  if(! length(b)) {
    xb <- 0 * xa
    db <- da
  }
  mb <- nrow(xb)

  if(length(a2)) {
    if(! length(b) || ! length(b2))
      stop('b and b2 must be given if a2 is given')
    da2 <- do.call('gendata', list(fit, factors=a2, expand=expand))
    xa2 <- pred(da2)
    ma2 <- nrow(xa2)
    db2 <- do.call('gendata', list(fit, factors=b2, expand=expand))
    xb2 <- pred(db2)
    mb2 <- nrow(xb2)
  }

  allsame <- function(x) diff(range(x)) == 0
  
  vary <- NULL
  mall <- c(ma, mb)
  ncols <- c(ncol(da), ncol(db))
  if(length(a2)) {
    mall <- c(mall, ma2, mb2)
    ncols <- c(ncols, ncol(da2), ncol(db2))
  }
  
  if(allsame(mall) && ! allsame(ncols)) stop('program logic error')
  if(any(sort(names(da)) != sort(names(db))))
    stop('program logic error')
  if(length(a2) && (any(sort(names(da)) != sort(names(da2))) ||
                       any(sort(names(da)) != sort(names(db2)))))
    stop('program logic error')
    
  if(TRUE) {
    ## If all lists have same length, label contrasts by any variable
    ## that has the same length and values in all lists
    k <- integer(0)
    nam <- names(da)
    for(j in 1 : length(da)) {
      w <- nam[j]
      eq <- all(as.character(da[[w]]) == as.character(db[[w]]))
      if(length(a2))
        eq <- eq & all(as.character(da[[w]]) == as.character(da2[[w]])) &
          all(as.character(db[[w]]) == as.character(db2[[w]]))
      if(eq) k <- c(k, j)
    }
    if(length(k)) vary <- da[k]
  } else if(max(mall) > 1) {
    ## Label contrasts by values of longest variable in list if
    ## it has the same length as the expanded design matrix
    d <- if(ma > 1) a else b
    if(length(a2) && (max(ma2, mb2) > max(ma, mb)))
      d <- if(ma2 > 1) a2 else b2
    l <- sapply(d, length)
    vary <- if(sum(l == max(mall)) == 1) d[l == max(mall)]
  }

  if(sum(mall > 1) > 1 && ! allsame(mall[mall > 1]))
    stop('lists of settings with more than one row must all have the same # rows')

  mm <- max(mall)
  if(mm > 1 && any(mall == 1)) {
    if(ma == 1) xa <- matrix(xa, nrow=mm, ncol=ncol(xa), byrow=TRUE)
    if(mb == 1) xb <- matrix(xb, nrow=mm, ncol=ncol(xb), byrow=TRUE)
    if(length(a2)) {
      if(ma2 == 1) xa2 <- matrix(xa2, nrow=mm, ncol=ncol(xa2), byrow=TRUE)
      if(mb2 == 1) xb2 <- matrix(xb2, nrow=mm, ncol=ncol(xb2), byrow=TRUE)
    }
  }
  
  X <- xa - xb
  if(length(a2)) X <- X - (xa2 - xb2)
  m <- nrow(X)
  
  if(is.character(weights)) {
    if(weights != 'equal') stop('weights must be "equal" or a numeric vector')
    weights <- rep(1,  m)
  } else if(length(weights) > 1)
      stop('can specify more than one weight only for unimplemented type="average"')
    else if(length(weights) != m) stop(paste('there must be', m, 'weights'))
  weights <- as.vector(weights)
  if(m > 1)
    X <- matrix(apply(weights*X, 2, sum) / sum(weights), nrow=1,
                dimnames=list(NULL, dimnames(X)[[2]]))

  X
  }
