\name{residuals.ols}
\alias{residuals.ols}
\title{Residuals for ols}

\description{Computes various residuals and measures of influence for a
  fit from \code{ols}.}
\usage{
\method{residuals}{ols}(object, 
      type=c("ordinary", "score", "dfbeta", "dfbetas", 
             "dffit", "dffits", "hat", "hscore"), \dots)
}
\arguments{
\item{object}{
object created by \code{ols}.  Depending on \code{type}, you may have had to
specify \code{x=TRUE} to \code{ols}.
}
\item{type}{
type of residual desired.  \code{"ordinary"} refers to the usual residual.
\code{"score"} is the matrix of score residuals (contributions to first
derivative of log likelihood).
\code{dfbeta} and \code{dfbetas} mean respectively the raw and normalized matrix 
of changes in regression coefficients after
deleting in turn each observation.  The coefficients are normalized by their
standard errors.  \code{hat} contains the leverages --- diagonals of the ``hat'' matrix.
\code{dffit} and \code{dffits} contain respectively the difference and normalized
difference in predicted values when each observation is omitted. 
The S \code{lm.influence} function is used.  When \code{type="hscore"}, the
ordinary residuals are divided by one minus the corresponding hat
matrix diagonal element to make residuals have equal variance.
}
\item{\dots}{ignored}
}
\value{
a matrix or vector, with places for observations that were originally
deleted by \code{ols} held by \code{NA}s 
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{lm.influence}}, \code{\link{ols}},
\code{\link{which.influence}} 
}
\examples{
set.seed(1)
x1 <- rnorm(100)
x2 <- rnorm(100)
x1[1] <- 100
y <- x1 + x2 + rnorm(100)
f <- ols(y ~ x1 + x2, x=TRUE, y=TRUE)
resid(f, "dfbetas")
which.influence(f)
}
\keyword{models}
\keyword{regression}
\concept{model validation}
