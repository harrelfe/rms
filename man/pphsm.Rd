\name{pphsm}
\alias{pphsm}
\alias{print.pphsm}
\alias{vcov.pphsm}
\title{Parametric Proportional Hazards form of AFT Models}

\description{
  Translates an accelerated failure time (AFT) model fitted by
  \code{psm} to proportional hazards form, if the fitted model was
  a Weibull or exponential model (extreme value distribution with
  "log" link).
}
\usage{
pphsm(fit)
\method{print}{pphsm}(x, digits=max(options()$digits - 4, 3),
correlation=TRUE, \dots)
\method{vcov}{pphsm}(object, \dots)
}
\arguments{
  \item{fit}{fit object created by \code{psm}}
  \item{x}{result of \code{psm}}
  \item{digits}{how many significant digits are to be used for the
    returned value}
  \item{correlation}{set to \code{FALSE} to suppress printing of
    correlation matrix of parameter estimates}
  \item{\dots}{ignored}
  \item{object}{a pphsm object}
}
\value{
a new fit object with transformed parameter estimates
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{psm}}, \code{\link{summary.rms}}, \code{\link{print.pphsm}}
}
\examples{
require(survival)
set.seed(1)
S <- Surv(runif(100))
x <- runif(100)
dd <- datadist(x); options(datadist='dd')
f <- psm(S ~ x, dist="exponential")
summary(f)        # effects on log(T) scale
f.ph <- pphsm(f)
\dontrun{summary(f.ph)     # effects on hazard ratio scale}
options(datadist=NULL)
}
\keyword{models}
\keyword{survival}
\keyword{regression}
