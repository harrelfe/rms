#Uses matinv Fortran function, which uses ginv and sweep
#Returns matrix inverse with attributes rank (integer rank of x)
# and swept (logical - whether or not ith variable has been swept)
#Input matrix should set its swept attribute before the first invocation
# of matinv for that matrix.  If swept isn't set, it defaults to all F.
#
#Inverse is with respect to diagonal elements which[1],which[2],...
#For collinearities, the appropriate rows and columns of inv are set to 0
#Caller must negate matrix when finished with all partial inversions if
# negate is false.  The default is to automatically negate the which
# portion of the inverse, i.e., to assume that no further operations are
# to be done on the matrix 
#
#Eps is singularity criterion, like 1-Rsquare
#
#F. Harrell 1 Aug 90

matinv <- function(a, which, negate=TRUE, eps=1E-12)
{
  swept <- attr(a,"swept")
  if(!is.matrix(a)) a <- as.matrix(a)
  storage.mode(a) <- "double"
  m<-nrow(a)
  if(missing(which))which <- 1:m
  else
    {
      rw <- range(which)
      if(rw[1] < 1 | rw[2] > m) stop("illegal elements to invert")
    }
  storage.mode(which) <- "integer"
  if(!length(swept))swept <- rep(FALSE, m)
  if(m!=ncol(a))stop("matrix must be square")

	y <- 
      .Fortran("matinv",x = a, as.integer(m), 
               as.integer(length(which)),which,
               swept=swept, logical(m), double(m*(m+1)/2), 
               double(m), rank = integer(1), as.double(eps),
               as.logical(negate), PACKAGE="rms")

  x <- y$x
  attr(x,"rank") <- y$rank
  attr(x,"swept") <- y$swept
  dimnames(x) <- dimnames(a)
  x
  
}
