\name{hazard.ratio.plot}
\alias{hazard.ratio.plot}
\title{Hazard Ratio Plot}

\description{
The \code{hazard.ratio.plot} function repeatedly estimates Cox
regression coefficients and confidence limits within time intervals.
The log hazard ratios are plotted against the mean failure/censoring
time within the interval. Unless \code{times} is specified, the number of
time intervals will be \eqn{\max(round(d/e),2)}, where \eqn{d} is the
total number 
of events in the sample. Efron's likelihood is used for estimating
Cox regression coefficients (using \code{coxph.fit}).  In the case of
tied failure times, some intervals may have a point in common.
}
\usage{
hazard.ratio.plot(x, Srv, which, times=, e=30, subset,
                  conf.int=.95, legendloc=NULL, smooth=TRUE, pr=FALSE, pl=TRUE,
                  add=FALSE, ylim, cex=.5, xlab="t", ylab, antilog=FALSE, \dots)
}
\arguments{
\item{x}{
a vector or matrix of predictors
}
\item{Srv}{a \code{Srv} or \code{Surv} object}
\item{which}{
a vector of column numbers of \code{x} for which to estimate hazard
ratios across time and make plots.
The default is to do so for all predictors.  Whenever
one predictor is displayed, all other predictors in the \code{x} matrix
are adjusted for (with a separate adjustment form for each time interval).
}
\item{times}{
optional vector of time interval endpoints.
Example: \code{times=c(1,2,3)} uses intervals \code{[0,1), [1,2), [2,3), [3+)}.
If times is omitted, uses intervals containing \code{e} events
}
\item{e}{
number of events per time interval if times not given
}
\item{subset}{
vector used for subsetting the entire analysis,
 e.g. \code{subset=sex=="female"}
}
\item{conf.int}{
confidence interval coverage
}
\item{legendloc}{
location for legend. Omit to use mouse, \code{"none"} for none,
 \code{"ll"} for lower left of graph, or actual x and y coordinates (e.g.
\code{c(2,3)})
}
\item{smooth}{
also plot the super--smoothed version of the log hazard ratios
}
\item{pr}{
defaults to \code{FALSE} to suppress printing of individual Cox fits
}
\item{pl}{
defaults to \code{TRUE} to plot results
}
\item{add}{
add this plot to an already existing plot
}
\item{ylim}{
vector of \code{y}-axis limits. Default is computed to include confidence bands.
}
\item{cex}{
character size for legend information, default is 0.5
}
\item{xlab}{
label for \code{x}-axis, default is \code{"t"}
}
\item{ylab}{
label for \code{y}-axis, default is \code{"Log Hazard Ratio"} or \code{"Hazard Ratio"},
depending on \code{antilog}.
}
\item{antilog}{
default is \code{FALSE}. Set to \code{TRUE} to plot anti-log, i.e., hazard ratio.
}
\item{...}{
optional graphical parameters
}}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link[survival]{cox.zph}}, \code{\link{residuals.cph}},
\code{\link[survival]{survival-internal}}, \code{\link{cph}},
\code{\link[survival]{coxph}}, \code{\link[survival]{Surv}}, \code{\link{Srv}}
}
\examples{
n <- 500
set.seed(1)
age <- 50 + 12*rnorm(n)
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50))
d.time <- -log(runif(n))/h
label(d.time) <- 'Follow-up Time'
e <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)
units(d.time) <- "Year"
hazard.ratio.plot(age, Srv(d.time,e), e=20, legendloc='ll')
}
\keyword{survival}
