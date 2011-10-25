\name{matinv}
\alias{matinv}
\title{
Total and Partial Matrix Inversion using Gauss-Jordan Sweep Operator
}
\description{
  This function inverts or partially inverts a matrix using pivoting
  (the sweep operator).  It is useful for sequential model-building.
  }
\usage{
matinv(a, which, negate=TRUE, eps=1e-12)
}
\arguments{
\item{a}{
square matrix to invert or partially invert.  May have been inverted or
partially inverted previously by matinv, in which case its "swept"
attribute is updated.  Will un-invert if already inverted.
}
\item{which}{
vector of column/row numbers in a to invert.  Default is all, for total
inverse.
}
\item{negate}{
So that the algorithm can keep track of which pivots have been swept
as well as roundoff errors, it actually returns the negative of the
inverse or partial inverse.  By default, these elements are negated to
give the usual expected result.  Set negate=FALSE if you will be passing
the result right back into matinv, otherwise, negate the submatrix
before sending back to matinv.
}
\item{eps}{
singularity criterion
}}
\value{
a square matrix, with attributes "rank" and "swept".
}
\references{
Clarke MRB (1982).  Algorithm AS 178: The Gauss-Jordan sweep operator
with detection of collinearity.  Appl Statist 31:166--9.

Ridout MS, Cobb JM (1986).  Algorithm AS R78 : A remark on algorithm AS 178:
The Gauss-Jordan
sweep operator with detection of collinearity.  Appl Statist 38:420--2.
}
\seealso{\code{\link{lrm}}, \code{\link{solve}}}

\examples{
a      <- diag(1:3)
a.inv1 <- matinv(a, 1, negate=FALSE)	     #Invert with respect to a[1,1]
a.inv1
a.inv  <- -matinv(a.inv1, 2:3, negate=FALSE) #Finish the job
a.inv
solve(a)
}
\keyword{array}
