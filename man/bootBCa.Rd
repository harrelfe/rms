\name{bootBCa}
\alias{bootBCa}
\title{BCa Bootstrap on Existing Bootstrap Replicates}
\description{
This functions constructs an object resembling one produced by the
\code{boot} package's \code{boot} function, and runs that package's
\code{boot.ci} function to compute BCa and percentile confidence limits.
}
\usage{bootBCa(estimate, estimates, n, seed, conf.int = 0.95)}
\arguments{
  \item{estimate}{original whole-sample estimate}
  \item{estimates}{vector of bootstrap estimates}
  \item{n}{original number of observations}
  \item{seed}{\code{.Random.seem} in effect before bootstrap estimates
	were run}
  \item{conf.int}{confidence level}
}
\value{a list containing two elements that are each 2-vectors, named
  \code{np} (percentile) and \code{bca}}
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
X <- cbind(1, x1=.4, x2=.6)
est <- X %*% coef(b)
ests <- as.vector(X %*% t(b$boot.Coef))
bootBCa(est, ests, n=100, seed=seed)
}}
\keyword{bootstrap}

