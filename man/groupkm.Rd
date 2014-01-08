\name{groupkm}
\alias{groupkm}
\title{Kaplan-Meier Estimates vs. a Continuous Variable}
\description{
Function to divide \code{x} (e.g. age, or predicted survival at time
\code{u} created by \code{survest}) into \code{g} quantile groups, get
Kaplan-Meier estimates at time \code{u} (a scaler), and to return a
matrix with columns \code{x}=mean \code{x} in quantile, \code{n}=number
of subjects, \code{events}=no. events, and \code{KM}=K-M survival at
time \code{u}, \code{std.err} = s.e. of -log K-M.  Confidence intervals
are based on -log S(t).  Instead of supplying \code{g}, the user can
supply the minimum number of subjects to have in the quantile group
(\code{m}, default=50).  If \code{cuts} is given
(e.g. \code{cuts=c(0,.1,.2,\dots,.9,.1)}), it overrides \code{m} and
\code{g}.  Calls Therneau's \code{survfitKM} in the \code{survival}
package to get Kaplan-Meiers estimates and standard errors.  }
\usage{
groupkm(x, Srv, m=50, g, cuts, u, 
        pl=FALSE, loglog=FALSE, conf.int=.95, xlab, ylab,
        lty=1, add=FALSE, cex.subtitle=.7, \dots)
}
\arguments{
\item{x}{variable to stratify}
\item{Srv}{
a \code{Surv} object - n x 2 matrix containing survival
time and event/censoring 
1/0 indicator.  Units of measurement come from the "units" attribute
of the survival time variable.  "Day" is the default.
}
\item{m}{desired minimum number of observations in a group}
\item{g}{number of quantile groups}
\item{cuts}{actual cuts in \code{x}, e.g. \code{c(0,1,2)} to use [0,1), [1,2].
}
\item{u}{time for which to estimate survival}
\item{pl}{TRUE to plot results}
\item{loglog}{
set to \code{TRUE} to plot \code{log(-log(survival))} instead of survival
}
\item{conf.int}{
defaults to \code{.95} for 0.95 confidence bars.  Set to \code{FALSE} to suppress bars.
}
\item{xlab}{
if \code{pl=TRUE}, is x-axis label.  Default is \code{label(x)} or name of calling argument
}
\item{ylab}{
if \code{pl=TRUE}, is y-axis label.  Default is constructed from \code{u} and time \code{units}
attribute.
}
\item{lty}{
line time for primary line connecting estimates
}
\item{add}{
set to \code{TRUE} if adding to an existing plot
}
\item{cex.subtitle}{
character size for subtitle. Default is \code{.7}.  Use \code{FALSE} to
suppress subtitle. 
}
\item{...}{plotting parameters to pass to the plot and errbar functions}
}
\value{
matrix with columns named \code{x} (mean predictor value in interval), \code{n} (sample size
in interval), \code{events} (number of events in interval), \code{KM} (Kaplan-Meier
estimate), \code{std.err} (standard error of -log \code{KM})
}
\seealso{
\code{\link[survival]{survfit}}, \code{\link[Hmisc]{errbar}},
\code{\link[Hmisc]{cut2}}, \code{\link[survival]{Surv}},
\code{\link[Hmisc]{units}}
}
\examples{
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50))
d.time <- -log(runif(n))/h
label(d.time) <- 'Follow-up Time'
e <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)
units(d.time) <- "Year"
groupkm(age, Surv(d.time, e), g=10, u=5, pl=TRUE)
#Plot 5-year K-M survival estimates and 0.95 confidence bars by 
#decile of age.  If omit g=10, will have >= 50 obs./group.
}
\keyword{survival}
\keyword{nonparametric}
\concept{grouping}
\concept{stratification}
\concept{aggregation}
