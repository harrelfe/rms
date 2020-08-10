\name{ie.setup}
\alias{ie.setup}
\title{Intervening Event Setup}

\description{
Creates several new variables which help set up a dataset for modeling
with \code{cph} or \code{coxph} when there is a single binary time-dependent
covariable which turns on at a given time, and stays on.  This is
typical when analyzing the impact of an intervening event.
\code{ie.setup} creates a \code{Surv} object using the start time, stop time
format.  It also creates a binary indicator for the intervening event,
and a variable called \code{subs} that is useful when \code{attach}-ing
a dataframe. 
\code{subs} has observation numbers duplicated for subjects having an
intervening event, so those subject's baseline covariables (that are
not time-dependent) can be duplicated correctly.
}
\usage{
ie.setup(failure.time, event, ie.time, break.ties=FALSE)
}
\arguments{
\item{failure.time}{
a numeric variable containing the event or censoring times for the
terminating event
}
\item{event}{
a binary (0/1) variable specifying whether observations had the
terminating event (event=1) or were censored (event=0)
}
\item{ie.time}{
intervening event times.  For subjects having no intervening events,
the corresponding values of ie.time must be NA.
}
\item{break.ties}{
Occasionally intervening events are recorded as happening at exactly
the same time as the termination of follow-up for some subjects.
The \code{Surv} and \code{Surv} functions will not allow this.  To
randomly break the ties 
by subtracting a random number from such tied intervening event times,
specify \code{break.ties=TRUE}.  The random number is uniform between zero and
the minimum difference between any two untied \code{failure.time}s.
}}
\value{
a list with components \code{S, ie.status, subs, reps}.  \code{S} is a
\code{Surv}
object containing start and stop times for intervals of observation, 
along with event indicators.  \code{ie.status} is one if the intervening
event has occurred at the start of the interval, zero otherwise.
\code{subs} is a vector of subscripts that can be used to replicate other
variables the same way \code{S} was replicated.  \code{reps} specifies how many
times each original observation was replicated.  \code{S, ie.status, subs} are
all the same length (at least the number of rows for \code{S} is) and are longer than the original \code{failure.time} vector.
\code{reps} is the same length as the original \code{failure.time} vector.
The \code{subs} vector is suitable for passing to \code{validate.lrm} or \code{calibrate},
which pass this vector under the name \code{cluster} on to \code{predab.resample} so that bootstrapping can be
done by sampling with replacement from the original subjects rather than
from the individual records created by \code{ie.setup}.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{cph}}, \code{\link[survival]{coxph}},
\code{\link[survival]{Surv}}, \code{\link{cr.setup}},
\code{\link{predab.resample}} 
}
\examples{
failure.time <- c(1 ,   2,   3)
event        <- c(1 ,   1,   0)
ie.time      <- c(NA, 1.5, 2.5)

z <- ie.setup(failure.time, event, ie.time)
S <- z$S
S
ie.status <- z$ie.status
ie.status
z$subs
z$reps
\dontrun{
attach(input.data.frame[z$subs,])   #replicates all variables
f <- cph(S ~ age + sex + ie.status)
# Instead of duplicating rows of data frame, could do this:
attach(input.data.frame)
z <- ie.setup(failure.time, event, ie.time)
s <- z$subs
age <- age[s]
sex <- sex[s]
f <- cph(S ~ age + sex + ie.status)
}
}
\keyword{survival}





