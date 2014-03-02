\name{vif}
\alias{vif}
\title{Variance Inflation Factors}

\description{
Computes variance inflation factors from the covariance matrix of
parameter estimates, using the method of Davis et al. (1986), which
is based on the correlation matrix from the information matrix.
}
\usage{
vif(fit)
}
\arguments{
\item{fit}{
an object created by \code{lrm}, \code{ols}, \code{psm}, \code{cph},
\code{Rq}, \code{Glm}, \code{glm} 
}}
\value{vector of vifs}
\author{
Frank Harrell
\cr
Department of Biostatistics
\cr
Vanderbilt University
\cr
f.harrell@vanderbilt.edu
}
\references{
Davis CE, Hyde JE, Bangdiwala SI, Nelson JJ: An example of dependencies 
among variables in a conditional logistic regression. In Modern
Statistical Methods in Chronic Disease Epidemiology, Eds SH Moolgavkar and
RL Prentice, pp. 140--147.  New York: Wiley; 1986.
}
\seealso{
\code{\link{rmsMisc}} (for \code{\link[Hmisc]{num.intercepts}}
}
\examples{
set.seed(1)
x1 <- rnorm(100)
x2 <- x1+.1*rnorm(100)
y  <- sample(0:1, 100, TRUE)
f  <- lrm(y ~ x1 + x2)
vif(f)
}
\keyword{models}
\keyword{regression}

