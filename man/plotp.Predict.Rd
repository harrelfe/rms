\name{plotp.Predict}
\alias{plotp.Predict}
\title{Plot Effects of Variables Estimated by a Regression Model Fit
	Using plotly}
\description{
  Uses \code{plotly} graphics (without using ggplot2) to plot the effect
  of one or two predictors 
  on the linear predictor or X beta scale, or on some transformation of
  that scale.  The first argument specifies the result of the
  \code{Predict} function.  The predictor is always plotted in its
  original coding.  Hover text shows point estimates, confidence
  intervals, and on the leftmost x-point, adjustment variable settings.

	If \code{Predict} was run with no variable settings, so that each
	predictor is varied one at a time, the result of \code{plotp.Predict}
	is a list with two elements.  The first, named \code{Continuous}, is a
	\code{plotly} object containing a single graphic with all the
	continuous predictors varying.  The second, named \code{Categorical},
	is a \code{plotly} object containing a single graphic with all the
	categorical predictors varying.  If there are no categorical
	predictors, the value returned by by \code{plotp.Predict} is a single
	\code{plotly} object and not a list of objects.

  If \code{rdata} is given, a spike histogram is drawn showing
  the location/density of data values for the \eqn{x}-axis variable.  If
  there is a superposition variable that generated separate
  curves, the data density specific to each class of points is shown.
  The histograms are drawn by \code{histSpikeg}.

  To plot effects instead of estimates (e.g., treatment differences as a
  function of interacting factors) see \code{contrast.rms} and
  \code{summary.rms}.

	Unlike \code{ggplot.Predict}, \code{plotp.Predict} does not handle
	\code{groups}, \code{anova}, or \code{perim} arguments.
}
\usage{
\method{plotp}{Predict}(data, subset, xlim, ylim, xlab, ylab, 
     rdata=NULL, nlevels=3, vnames=c('labels','names'),
     histSpike.opts=list(frac=function(f) 0.01 + 
         0.02 * sqrt(f - 1)/sqrt(max(f, 2) - 1), side=1, nint=100),
     ncols=3, width=800, ...)
}
\arguments{
\item{data}{a data frame created by \code{Predict}}
\item{subset}{a subsetting expression for restricting the rows of
  \code{data} that are used in plotting.  For example, predictions may have
  been requested for males and females but one wants to plot only females.}
\item{xlim}{ignored unless predictors were specified to \code{Predict}.
	Specifies the x-axis limits of the single plot produced.}
\item{ylim}{
Range for plotting on response variable axis. Computed by default and
includes the confidence limits.
}
\item{xlab}{
Label for \code{x}-axis when a single plot is made, i.e., when a
predictor is specified to \code{Predict}.  Default is one given to
\code{asis, rcs}, etc., which may have been the \code{"label"} attribute
of the variable. 
}
\item{ylab}{
Label for \code{y}-axis.  If \code{fun} is not given,
default is \code{"log Odds"} for
\code{lrm}, \code{"log Relative Hazard"} for \code{cph}, name of the response
variable for \code{ols}, \code{TRUE} or \code{log(TRUE)} for \code{psm},
  or \code{"X * Beta"} otherwise.  Specify \code{ylab=NULL} to omit
  \code{y}-axis labels.
}
\item{rdata}{a data frame containing the original raw data on which the
  regression model were based, or at least containing the \eqn{x}-axis
  and grouping variable.  If \code{rdata} is present and contains the
  needed variables, the original data are added to the graph in the form
  of a spike histogram using \code{histSpikeg} in the Hmisc package.
}
\item{nlevels}{
  A non-numeric x-axis variable with \code{nlevels} or fewer unique values
  will cause a horizontal dot plot to be drawn instead of an x-y plot.
}
\item{vnames}{applies to the case where multiple plots are produced
  separately by predictor.  Set to \code{'names'} to use variable names
  instead of labels for these small plots.}
\item{histSpike.opts}{a list containing named elements that specifies
  parameters to \code{\link[Hmisc:scat1d]{histSpikeg}} when \code{rdata} is given.  The
  \code{col} parameter is usually derived from other plotting
  information and not specified by the user.}
\item{ncols}{number of columns of plots to use when plotting multiple
	continuous predictors}
\item{width}{width in pixels for \code{plotly} graphics}
\item{\dots}{ignored}
}
\value{a \code{plotly} object or a list containing two elements, each
	one a \code{plotly} object}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\references{
  Fox J, Hong J (2009): Effect displays in R for multinomial and
  proportional-odds logit models: Extensions to the effects package.  J
  Stat Software 32 No. 1.
}
\seealso{
  \code{\link{Predict}}, \code{\link{rbind.Predict}},
  \code{\link{datadist}}, \code{\link{predictrms}},
  \code{\link{contrast.rms}}, \code{\link{summary.rms}},
  \code{\link{rms}}, \code{\link{rmsMisc}}, \code{\link{plot.Predict}},
	\code{\link{ggplot.Predict}},
  \code{\link[Hmisc]{histSpikeg}},
  \code{\link[Hmisc]{Overview}}
}
\examples{
\dontrun{
n <- 350     # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
    (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male')) +
    .01 * (blood.pressure - 120)
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
           x=TRUE, y=TRUE)

p <- plotp(Predict(fit))
p$Continuous
p$Categorical
# When using Rmarkdown html notebook, best to use
# prList(p) to render the two objects
plotp(Predict(fit), rdata=llist(blood.pressure, age))$Continuous
# spike histogram plot for two of the predictors

p <- Predict(fit, name=c('age','cholesterol'))   # Make 2 plots
plotp(p)

p <- Predict(fit, age, sex)
plotp(p, rdata=llist(age,sex))
# rdata= allows rug plots (1-dimensional scatterplots)
# on each sex's curve, with sex-
# specific density of age
# If data were in data frame could have used that
p <- Predict(fit, age=seq(20,80,length=100), sex='male', fun=plogis)
# works if datadist not used
plotp(p, ylab='P')
# plot predicted probability in place of log odds

# Compute predictions for three predictors, with superpositioning or
# conditioning on sex, combined into one graph

p1 <- Predict(fit, age, sex)
p2 <- Predict(fit, cholesterol, sex)
p3 <- Predict(fit, blood.pressure, sex)
p <- rbind(age=p1, cholesterol=p2, blood.pressure=p3)
plotp(p, ncols=2, rdata=llist(age, cholesterol, sex))
}
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
