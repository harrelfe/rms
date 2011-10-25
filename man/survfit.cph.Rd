\name{survfit.cph}
\alias{survfit.cph}
\title{
Cox Predicted Survival
}
\description{
This is a slightly modified version of Therneau's \code{survfit.coxph}
function.  The difference is that \code{survfit.cph} assumes that
\code{x=TRUE,y=TRUE} were specified to the fit. This assures that the
environment in effect at the time of the fit (e.g., automatic knot
estimation for spline functions) is the same one used for basing predictions.
}
\usage{
\method{survfit}{cph}(formula, newdata, se.fit=TRUE, conf.int=0.95, 
        individual=FALSE, type=NULL, vartype=NULL,
        conf.type=c('log', "log-log", "plain", "none"), id, \dots)
}
\arguments{
\item{formula}{
a fit object from \code{cph} or \code{coxph}
see \code{\link[survival]{survfit.coxph}}
}
\item{newdata,se.fit,conf.int,individual,type,vartype,conf.type,id}{see
  \code{\link[survival]{survfit}}.  If \code{individual} is \code{TRUE},
  there must be exactly one \code{Surv} object in \code{newdata}.  This
  object is used to specify time intervals for time-dependent covariate
  paths.  To get predictions for multiple subjects with time-dependent
  covariates, specify a vector \code{id} which specifies unique
  hypothetical subjects.  The length of \code{id} should equal the
  number of rows in \code{newdata}.}
\item{\dots}{Not used}
}
\value{see \code{survfit.coxph}}
\seealso{\code{\link{survest.cph}}}
\keyword{survival}
