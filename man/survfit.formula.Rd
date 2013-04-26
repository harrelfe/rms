\name{survfit.formula}
\alias{survfit.formula}
\title{Compute a Survival Curve for Censored Data}
\description{
Computes an estimate of a survival curve for censored data
using either the Kaplan-Meier or the Fleming-Harrington method
or computes the predicted survivor function.
For competing risks data it computes the cumulative incidence curve.
This calls the \code{survival} package's \code{survfit.formula} function
with a different default value for \code{conf.type} (log-log basis).  In
addition, attributes of the event time variable are saved (label and
units of measurement).
}
\usage{
\method{survfit}{formula}(formula, data, \dots)
}
\arguments{
  \item{formula}{
    a formula object, which must have a \code{Srv} or \code{Surv} object as the
    response on the left of the \code{~} operator and, if desired, terms
    separated by + operators on the right.
    One of the terms may be a \code{strata} object.
    For a single survival curve the right hand side should be \code{~ 1}.
  }
  \item{data}{
    a data frame in which to interpret the variables named in the formula,
    \code{subset} and \code{weights} arguments.
  }
  \item{\dots}{see \code{\link[survival]{survfit.formula}}}
}
\value{
an object of class \code{"survfit"}.
See \code{survfit.object} for
details. Methods defined for survfit objects are
\code{print}, \code{plot},
\code{lines}, and \code{points}.
}
\details{
  see \code{\link[survival]{survfit.formula}} for details
}
\seealso{
\code{\link{survfit.cph}} for survival curves from Cox models.
\code{\link{print}},
\code{\link{plot}},
\code{\link{lines}},
\code{\link{coxph}},
\code{\link{Srv}},
\code{\link{strata}}.
}
\author{Thomas Lumley \email{tlumley@u.washington.edu} and Terry Therneau}

\examples{
require(survival)
#fit a Kaplan-Meier and plot it
fit <- survfit(Srv(time, status) ~ x, data = aml)
plot(fit, lty = 2:3)
legend(100, .8, c("Maintained", "Nonmaintained"), lty = 2:3)

#fit a Cox proportional hazards model and plot the
#predicted survival for a 60 year old
fit <- coxph(Srv(futime, fustat) ~ age, data = ovarian)
plot(survfit(fit, newdata=data.frame(age=60)),
     xscale=365.25, xlab = "Years", ylab="Survival")

# Here is the data set from Turnbull
#  There are no interval censored subjects, only left-censored (status=3),
#  right-censored (status 0) and observed events (status 1)
#
#                             Time
#                         1    2   3   4
# Type of observation
#           death        12    6   2   3
#          losses         3    2   0   3
#      late entry         2    4   2   5
#
tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
                    status=rep(c(1,0,2),4),
                    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))
fit  <- survfit(Srv(time, time, status, type='interval') ~1,
                data=tdata, weights=n)

#
# Time to progression/death for patients with monoclonal gammopathy
#  Competing risk curves (cumulative incidence)
fit1 <- survfit(Srv(stop, event=='progression') ~1, data=mgus1,
                    subset=(start==0))
fit2 <- survfit(Srv(stop, status) ~1, data=mgus1,
                    subset=(start==0), etype=event) #competing risks
# CI curves are always plotted from 0 upwards, rather than 1 down
plot(fit2, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            col=2:3, xlab="Years post diagnosis of MGUS")
lines(fit1, fun='event', xscale=365.25, xmax=7300, mark.time=FALSE,
            conf.int=FALSE)
text(10, .4, "Competing Risk: death", col=3)
text(16, .15,"Competing Risk: progression", col=2)
text(15, .30,"KM:prog")
}
