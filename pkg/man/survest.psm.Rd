\name{survest.psm}
\alias{survest.psm}
\alias{print.survest.psm}
\title{Parametric Survival Estimates}
\description{
Computes predicted survival probabilities or hazards and optionally confidence
limits (for survival only) for parametric survival models fitted with
\code{psm}. 
If getting predictions for more than one observation, \code{times} must
be specified. For a model without predictors, no input data are
specified.
}
\usage{
\method{survest}{psm}(fit, newdata, linear.predictors, x, times, fun,
        loglog=FALSE, conf.int=0.95,
        what=c("survival","hazard","parallel"), \dots)

\method{print}{survest.psm}(x, \dots)
}
\arguments{
\item{fit}{
fit from \code{psm}
}
\item{newdata, linear.predictors, x, times, conf.int}{
see \code{survest.cph}. One of \code{newdata}, \code{linear.predictors}, \code{x} must be given.
\code{linear.predictors} includes the intercept.
If \code{times} is omitted, predictions are made at 200 equally spaced points
between 0 and the maximum failure/censoring time used to fit the model.

\code{x} can also be a result from \code{survest.psm}.
}
\item{what}{
The default is to compute survival probabilities.  Set \code{what="hazard"} or
some abbreviation of \code{"hazard"} to compute hazard rates.
\code{what="parallel"} assumes that the length of \code{times} is the number of
subjects (or one), and causes \code{survest} to estimate the
\eqn{i^{th}} subject's survival probability at the \eqn{i^{th}} value of
\code{times} (or at the scalar value of \code{times}). 
\code{what="parallel"} is used by \code{val.surv} for example.
}
\item{loglog}{
set to \code{TRUE} to transform survival estimates and confidence limits using
log-log
}
\item{fun}{
a function to transform estimates and optional confidence intervals
}
\item{\dots}{unused}
}
\value{
see \code{survest.cph}. If the model has no predictors, predictions are
made with respect to varying time only, and the returned object
is of class \code{"survfit"} so the survival curve can be plotted
with survplot.survfit. If \code{times} is omitted, the
entire survival curve or hazard from \code{t=0,\dots,fit$maxtime} is estimated, with
increments computed to yield 200 points where \code{fit$maxtime} is the
maximum survival time in the data used in model fitting. Otherwise,
the \code{times} vector controls the time points used.
}
\details{
Confidence intervals are based on asymptotic normality of the linear
predictors.
The intervals account for the fact that a scale parameter may have been
estimated jointly with beta.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{psm}}, \code{\link[survival]{survreg}}, \code{\link{rms}}, \code{\link[survival]{survfit}}, \code{\link{predictrms}}, \code{\link{survplot}},
\code{\link[survival]{survreg.distributions}}
}
\examples{
# Simulate data from a proportional hazards population model
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
S <- Surv(dt,e)

f <- psm(S ~ lsp(age,c(40,70)))
survest(f, data.frame(age=seq(20,80,by=5)), times=2)

#Get predicted survival curve for 40 year old
survest(f, data.frame(age=40))

#Get hazard function for 40 year old
survest(f, data.frame(age=40), what="hazard")$surv #still called surv
}
\keyword{survival}
\keyword{regression}
\keyword{models}
