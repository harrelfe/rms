\name{survest.cph}
\alias{survest}
\alias{survest.cph}
\title{
Cox Survival Estimates
}
\description{
Compute survival probabilities and optional confidence limits for
Cox survival models.  If \code{x=TRUE, y=TRUE} were specified to \code{cph},
confidence limits use the correct formula for any combination of
predictors. Otherwise, if \code{surv=TRUE} was specified to \code{cph},
confidence limits are based only on standard errors of \code{log(S(t))}
at the mean value of \eqn{X\beta}{X beta}. If the model 
contained only stratification factors, or if predictions are being
requested near the mean of each covariable, this approximation will be
accurate. Unless \code{times} is given, at most one observation may be
predicted.
}
\usage{
survest(fit, \dots)
\method{survest}{cph}(fit, newdata, linear.predictors, x, times, 
        fun, loglog=FALSE, conf.int=0.95, type, vartype,
        conf.type=c("log", "log-log", "plain", "none"), se.fit=TRUE,
        what=c('survival','parallel'),
        individual=FALSE, ...)
}
\arguments{
\item{fit}{
a model fit from \code{cph}
}
\item{newdata}{
a data frame containing predictor variable combinations for which
predictions are desired
}
\item{linear.predictors}{
a vector of linear predictor values (centered) for which predictions
are desired. If the model is stratified, the "strata" attribute
must be attached to this vector (see example).
}
\item{x}{
 a design matrix at which to compute estimates, with any strata attached
 as a "strata" attribute. Only one of \code{newdata},
 \code{linear.predictors}, or \code{x} may be specified.  If none is
 specified, but \code{times} is specified, you will get survival
 predictions at all subjects' linear predictor and strata values.
}
\item{times}{
a vector of times at which to get predictions. If omitted, predictions
are made at all unique failure times in the original input data.
}
\item{loglog}{
set to \code{TRUE} to make the \code{log-log} transformation of survival
estimates and confidence limits.
}
\item{fun}{
any function to transform the estimates and confidence limits (\code{loglog}
is a special case)
}
\item{conf.int}{
set to \code{FALSE} or \code{0} to suppress confidence limits, or e.g. \code{.95} to 
cause 0.95 confidence limits to be computed
}
\item{type}{
see \code{survfit.coxph}
}
\item{vartype}{
see \code{survfit.coxph}
}
\item{conf.type}{
specifies the basis for computing confidence limits. \code{"log"} is the
default as in the \code{survival} package.
}
\item{se.fit}{
set to \code{TRUE} to get standard errors of log predicted survival
(no matter what \code{conf.type} is).
If \code{FALSE}, confidence limits are suppressed.
}
\item{individual}{
set to \code{TRUE} to have \code{survfit} interpret \code{newdata} as
specifying a covariable path for a single individual (represented by
multiple records).
}
\item{what}{
Normally use \code{what="survival"} to estimate survival probabilities at
times that may not correspond to the subjects' own times.
\code{what="parallel"} assumes that the length of \code{times} is the number of
subjects (or one), and causes \code{survest} to estimate the ith subject's survival
probability at the ith value of \code{times} (or at the scalar value of \code{times}).
\code{what="parallel"} is used by \code{val.surv} for example.
}
\item{\dots}{unused}
}
\value{
  If \code{times} is omitted, returns a list with the elements
  \code{time}, \code{n.risk}, \code{n.event}, \code{surv}, \code{call}
  (calling statement), and optionally \code{std.err}, \code{upper},
  \code{lower}, \code{conf.type}, \code{conf.int}. The estimates in this
  case correspond to one subject. If \code{times} is specified, the
  returned list has possible components \code{time}, \code{surv},
  \code{std.err}, \code{lower}, and \code{upper}. These will be matrices
  (except for \code{time}) if more than one subject is being predicted,
  with rows representing subjects and columns representing \code{times}.
  If \code{times} has only one time, these are reduced to vectors with
  the number of elements equal to the number of subjects.  }
\details{
  The result is passed through \code{naresid} if \code{newdata},
  \code{linear.predictors}, and \code{x} are not specified, to restore
  placeholders for \code{NA}s.
  }
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{cph}}, \code{\link{survfit.cph}}, \code{\link[survival]{survfit.coxph}}, \code{\link{predictrms}}, \code{\link{survplot}}
}
\examples{
# Simulate data from a population model in which the log hazard
# function is linear in age and there is no age x sex interaction
# Proportional hazards holds for both variables but we
# unnecessarily stratify on sex to see what happens
require(survival)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('Male','Female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
Srv <- Surv(dt,e)


f <- cph(Srv ~ age*strat(sex), x=TRUE, y=TRUE) #or surv=T
survest(f, expand.grid(age=c(20,40,60),sex=c("Male","Female")),
	    times=c(2,4,6), conf.int=.9)
f <- update(f, surv=TRUE)
lp <- c(0, .5, 1)
f$strata   # check strata names
attr(lp,'strata') <- rep(1,3)  # or rep('sex=Female',3)
survest(f, linear.predictors=lp, times=c(2,4,6))

# Test survest by comparing to survfit.coxph for a more complex model
f <- cph(Srv ~ pol(age,2)*strat(sex), x=TRUE, y=TRUE)
survest(f, data.frame(age=median(age), sex=levels(sex)), times=6)

age2 <- age^2
f2 <- coxph(Srv ~ (age + age2)*strata(sex))
new <- data.frame(age=median(age), age2=median(age)^2, sex='Male')
summary(survfit(f2, new), times=6)
new$sex <- 'Female'
summary(survfit(f2, new), times=6)

options(datadist=NULL)
}
\keyword{models}
\keyword{survival}
\keyword{regression}
