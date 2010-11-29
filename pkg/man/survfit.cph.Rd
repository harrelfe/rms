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
        individual=FALSE, type, vartype,
        conf.type=c('log', "log-log", "plain", "none"), \dots)
}
\arguments{
\item{formula}{
a fit object from \code{cph} or \code{coxph}
see \code{\link[survival]{survfit.coxph}}
}
\item{newdata,se.fit,conf.int,individual,type,vartype,conf.type}{see \code{\link[survival]{survfit}}}
\item{\dots}{Not used}
}
\value{see \code{survfit.coxph}}
\seealso{\code{\link{survest.cph}}}
\keyword{survival}
