\name{bootBCa}
\alias{bootBCa}
\title{BCa Bootstrap on Existing Bootstrap Replicates}
\description{
This functions constructs an object resembling one produced by the
\code{boot} package's \code{boot} function, and runs that package's
\code{boot.ci} function to compute BCa and percentile confidence limits.
\code{bootBCa} can provide separate confidence limits for a vector of
statistics when \code{estimate} has length greater than 1.  In that
case, \code{estimates} must have the same number of columns as
\code{estimate} has values.
}
\usage{bootBCa(estimate, estimates, type=c('percentile','bca','basic'),
               n, seed, conf.int = 0.95)}
\arguments{
  \item{estimate}{original whole-sample estimate}
  \item{estimates}{vector of bootstrap estimates}
  \item{type}{type of confidence interval, defaulting to nonparametric
	percentile}
  \item{n}{original number of observations}
  \item{seed}{\code{.Random.seem} in effect before bootstrap estimates
	were run}
  \item{conf.int}{confidence level}
}
\value{a 2-vector if \code{estimate} is of length 1, otherwise a matrix
  with 2 rows and number of columns equal to the length of
  \code{estimate}}
\author{Frank Harrell}
\note{
You can use \code{if(!exists('.Random.seed')) runif(1)} before running
your bootstrap to make sure that \code{.Random.seed} will be available
to \code{bootBCa}.
}
\seealso{\code{\link[boot]{boot.ci}}}
\examples{
\dontrun{
x1 <- runif(100); x2 <- runif(100); y <- sample(0:1, 100, TRUE)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
seed <- .Random.seed
b <- bootcov(f)
# Get estimated log odds at x1=.4, x2=.6
X <- cbind(c(1,1), x1=c(.4,2), x2=c(.6,3))
est <- X %*% coef(b)
ests <- t(X %*% t(b$boot.Coef))
bootBCa(est, ests, n=100, seed=seed)
bootBCa(est, ests, type='bca', n=100, seed=seed)
bootBCa(est, ests, type='basic', n=100, seed=seed)
}}
\keyword{bootstrap}

