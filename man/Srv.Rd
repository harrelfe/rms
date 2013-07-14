\name{Srv}
\alias{Srv}
\alias{[.Srv}
\title{Front-end to \code{survival} package \code{Surv} function}
\description{
  \code{Srv} saves labeling and time units information when creating
  \code{Surv} objects.  The subsetting method preserves these extra attributes.
}
\usage{
Srv(time, time2, event,
    type = c("right", "left", "interval", "counting", "interval2",
             "mstate"), origin = 0)
}
\arguments{
  \item{time,time2,event,type,origin}{see \code{\link[survival]{Surv}}}}
\value{a matrix of class \code{"Srv"} plus classes from \code{Surv}}
\author{Frank Harrell}
\seealso{\code{\link[survival]{Surv}}}
\examples{
dtime <- c(1, 3, 7); death <- c(1, 0, 1)
units(dtime) <- 'Month'
label(dtime) <- 'Event/censoring time'
label(death) <- 'Death, all causes'
S <- Srv(dtime, death)
attributes(S)
S[2:3,]
attributes(S[2:3,])
}
\keyword{survival}
