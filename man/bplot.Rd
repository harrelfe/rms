\name{bplot}
\alias{bplot}
\alias{perimeter}
\title{
3-D Plots Showing Effects of Two Continuous Predictors in a Regression
Model Fit}
\description{
Uses lattice graphics and the output from \code{Predict} to plot image,
contour, or perspective plots showing the simultaneous effects of two
continuous predictor variables.  Unless \code{formula} is provided, the
\eqn{x}-axis is constructed from the first variable listed in the call
to \code{Predict} and the \eqn{y}-axis variable comes from the second.

The \code{perimeter} function is used to generate the boundary of data
to plot when a 3-d plot is made.  It finds the area where there are
sufficient data to generate believable interaction fits.
}
\usage{
bplot(x, formula, lfun=lattice::levelplot, xlab, ylab, zlab,
      adj.subtitle=!info$ref.zero, cex.adj=.75, cex.lab=1,
      perim, showperim=FALSE,
      zlim=range(yhat, na.rm=TRUE), scales=list(arrows=FALSE),
      xlabrot, ylabrot, zlabrot=90, \dots)

perimeter(x, y, xinc=diff(range(x))/10, n=10, lowess.=TRUE)
}
\arguments{
  \item{x}{
	for \code{bplot}, an object created by \code{Predict} for which
    two or more numeric predictors varied.
	For \code{perim} is
	the first variable of a pair of predictors forming a 3-d plot.
  }
  \item{formula}{
	a formula of the form \code{f(yhat) ~ x*y} optionally followed by |a*b*c
	which are 1-3 paneling variables that were specified to \code{Predict}.
	\code{f} can represent any R function of a vector that produces a
	vector.  If the left hand side of the formula is omitted,
	\code{yhat} will be inserted.  If \code{formula} is omitted, it will
	be inferred from the first two variables that varied in the call to
	\code{Predict}.
  }
  \item{lfun}{
	a high-level lattice plotting function that takes formulas of the
	form \code{z ~ x*y}.  The default is an image plot
	(\code{levelplot}).  Other common choices are \code{wireframe} for
	perspective plot or \code{contourplot} for a contour plot.
	}
\item{xlab}{
  Character string label for \eqn{x}-axis. Default is given by \code{Predict}.
}
\item{ylab}{
  Character string abel for \eqn{y}-axis
}
\item{zlab}{
  Character string \eqn{z}-axis label for perspective (wireframe) plots.
  Default comes from \code{Predict}.  \code{zlab} will often be
  specified if \code{fun} was specified to \code{Predict}.
  }
\item{adj.subtitle}{
Set to \code{FALSE} to suppress subtitling the graph with the list of
       settings of non-graphed adjustment values. Default is \code{TRUE}
       if there are non-plotted adjustment variables and \code{ref.zero}
	   was not used.
}
\item{cex.adj}{
  \code{cex} parameter for size of adjustment settings in subtitles.  Default is
  0.75
}
\item{cex.lab}{
  \code{cex} parameter for axis labels.  Default is 1.
}
\item{perim}{
  names a matrix created by \code{perimeter} when used for 3-d plots of
  two continuous predictors.  When the combination of variables is outside
  the range in \code{perim}, that section of the plot is suppressed.  If
  \code{perim} 
  is omitted, 3-d plotting will use the marginal distributions of the
  two predictors to determine the plotting region, when the grid is
  not specified explicitly in \code{variables}.  When instead a series of
  curves is being plotted, \code{perim} specifies a function having two
  arguments.  The first is the vector of values of the first variable that
  is about to be plotted on the \eqn{x}-axis.  The second argument
  is the single 
  value of the variable representing different curves, for the current
  curve being plotted.  The function's returned value must be a logical
  vector whose length is the same as that of the first argument, with
  values \code{TRUE} if the corresponding point should be plotted for the
  current curve, \code{FALSE} otherwise.  See one of the latter examples.
}
\item{showperim}{
  set to \code{TRUE} if \code{perim} is specified and you want to
  show the actual perimeter used.
}
\item{zlim}{
  Controls the range for plotting in the \eqn{z}-axis if there is
  one. Computed by default.
}
\item{scales}{see \code{\link[lattice:cloud]{wireframe}}
}
\item{xlabrot}{rotation angle for the x-axis.  Default is 30 for
  \code{wireframe} and 0 otherwise.
  }
\item{ylabrot}{rotation angle for the y-axis.  Default is -40 for
  \code{wireframe}, 90 for \code{contourplot} or \code{levelplot},
  and 0 otherwise.
}
\item{zlabrot}{rotation angle for z-axis rotation for
  \code{wireframe} plots
}
\item{\dots}{other arguments to pass to the lattice function
}
\item{y}{
  second variable of the pair for \code{perim}.  If omitted, \code{x} is
  assumed to be a list with both \code{x} and \code{y} components.
}
\item{xinc}{
  increment in \code{x} over which to examine the density of \code{y} in
  \code{perimeter} 
}
\item{n}{
  within intervals of \code{x} for \code{perimeter}, takes the
  informative range of \code{y} to be the \eqn{n}th smallest to the
  \eqn{n}th largest values of \code{y}.  If there aren't 
  at least 2\eqn{n} \code{y} values in the \code{x} interval, no
  \code{y} ranges are used for that interval.
}
\item{lowess.}{
  set to \code{FALSE} to not have \code{lowess} smooth the data perimeters
}
}
\value{
\code{perimeter} returns a matrix of class \code{perimeter}.  This
       outline can be conveniently plotted by \code{lines.perimeter}.
}
\details{
\code{perimeter} is a kind of generalization of \code{datadist} for 2
continuous variables.  First, the \code{n} smallest and largest \code{x}
values are determined.  These form the lowest and highest possible
\code{x}s to display.  Then \code{x} is grouped into intervals bounded
by these two numbers, with the interval widths defined by \code{xinc}.
Within each interval, \code{y} is sorted and the \eqn{n}th smallest and
largest \code{y} are taken as the interval containing sufficient data
density to plot interaction surfaces.  The interval is ignored when
there are insufficient \code{y} values.  When the data are being
readied for \code{persp}, \code{bplot} uses the \code{approx} function to do
linear interpolation of the \code{y}-boundaries as a function of the
\code{x} values actually used in forming the grid (the values of the
first variable specified to \code{Predict}).  To make the perimeter smooth,
specify \code{lowess.=TRUE} to \code{perimeter}.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{datadist}}, \code{\link{Predict}},
\code{\link{rms}}, \code{\link{rmsMisc}}, \code{\link[lattice]{levelplot}},
       \code{\link[lattice]{contourplot}}, \code{\link[lattice:cloud]{wireframe}}
}
\examples{
n <- 1000    # define sample size
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
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
               x=TRUE, y=TRUE)
p <- Predict(fit, age, cholesterol, sex, np=50) # vary sex last
bplot(p)                 # image plot for age, cholesterol with color
                         # coming from yhat; use default ranges for
                         # both continuous predictors; two panels (for sex)
bplot(p, lfun=wireframe) # same as bplot(p,,wireframe)
# View from different angle, change y label orientation accordingly
# Default is z=40, x=-60
bplot(p,, wireframe, screen=list(z=40, x=-75), ylabrot=-25)
bplot(p,, contourplot)   # contour plot
bounds  <- perimeter(age, cholesterol, lowess=TRUE)
plot(age, cholesterol)     # show bivariate data density and perimeter
lines(bounds[,c('x','ymin')]); lines(bounds[,c('x','ymax')])
p <- Predict(fit, age, cholesterol)  # use only one sex
bplot(p, perim=bounds)   # draws image() plot
                         # don't show estimates where data are sparse
                         # doesn't make sense here since vars don't interact
bplot(p, plogis(yhat) ~ age*cholesterol) # Probability scale
options(datadist=NULL)
}
\keyword{models}
\keyword{hplot}
\keyword{htest}
