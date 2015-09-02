\name{survplot}
\alias{survplot}
\alias{survplot.rms}
\alias{survplot.npsurv}
\alias{survdiffplot}
\title{
Plot Survival Curves and Hazard Functions
}
\description{
Plot estimated survival curves, and for parametric survival models, plot
hazard functions.  There is an option to print the number of subjects
at risk at the start of each time interval.  Curves are automatically
labeled at the points of maximum separation (using the \code{labcurve}
function), and there are many other options for labeling that can be
specified with the \code{label.curves} parameter.  For example, different
plotting symbols can be placed at constant x-increments and a legend
linking the symbols with category labels can automatically positioned on
the most empty portion of the plot.

For the case of a two stratum analysis by \code{npsurv},
\code{survdiffplot} plots the difference in two Kaplan-Meier estimates
along with approximate confidence bands for the differences, with a
reference line at zero.  The number of subjects at risk is optionally
plotted.  This number is taken as the minimum of the number of subjects
at risk over the two strata.  When \code{conf='diffbands'},
\code{survdiffplot} instead does not make a new plot but adds a shaded
polygon to an existing plot, showing the midpoint of two survival
estimates plus or minus 1/2 the width of the confidence interval for the
difference of two Kaplan-Meier estimates.
}
\usage{
survplot(fit, \dots)
\method{survplot}{rms}(fit, \dots, xlim,
         ylim=if(loglog) c(-5, 1.5) else if
                 (what == "survival" & missing(fun)) c(0, 1),
         xlab, ylab, time.inc,
         what=c("survival","hazard"),
         type=c("tsiatis","kaplan-meier"),
         conf.type=c("log","log-log","plain","none"),
         conf.int=FALSE, conf=c("bands","bars"),
         add=FALSE, label.curves=TRUE,
         abbrev.label=FALSE, levels.only=FALSE,
         lty, lwd=par("lwd"),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         adj.subtitle=TRUE, loglog=FALSE, fun,
         n.risk=FALSE, logt=FALSE, dots=FALSE, dotsize=.003,
         grid=NULL, srt.n.risk=0, sep.n.risk=0.056, adj.n.risk=1, 
         y.n.risk, cex.n.risk=.6, pr=FALSE)
\method{survplot}{npsurv}(fit, xlim, 
         ylim, xlab, ylab, time.inc,
         conf=c("bands","bars","diffbands","none"), add=FALSE, 
         label.curves=TRUE, abbrev.label=FALSE,
         levels.only=FALSE, lty,lwd=par('lwd'),
         col=1, col.fill=gray(seq(.95, .75, length=5)),
         loglog=FALSE, fun, n.risk=FALSE, aehaz=FALSE, times=NULL,
         logt=FALSE, dots=FALSE, dotsize=.003, grid=NULL,
         srt.n.risk=0, sep.n.risk=.056, adj.n.risk=1,
         y.n.risk, cex.n.risk=.6, pr=FALSE, \dots)
survdiffplot(fit, order=1:2, fun=function(y) y,
           xlim, ylim, xlab, ylab="Difference in Survival Probability",
           time.inc, conf.int, conf=c("shaded", "bands","diffbands","none"),
           add=FALSE, lty=1, lwd=par('lwd'), col=1,
           n.risk=FALSE, grid=NULL,
           srt.n.risk=0, adj.n.risk=1,
           y.n.risk, cex.n.risk=.6)
}
\arguments{
\item{fit}{
result of fit (\code{cph}, \code{psm}, \code{npsurv},
\code{survest.psm}).  For \code{survdiffplot}, \code{fit} must be the
result of \code{npsurv}.
}
\item{\dots}{
list of factors with names used in model. For fits from \code{npsurv}
these arguments do not appear - all strata are plotted. Otherwise the
first factor listed is the factor used to determine different survival
curves.  Any other factors are used to specify single constants to be
adjusted to, when defaults given to fitting routine (through
\code{limits}) are not used.  The value given to factors is the original
coding of data given to fit, except that for categorical or strata
factors the text string levels may be specified.  The form of values
given to the first factor are none (omit the equal sign to use default range
or list of all values if variable is discrete), \code{"text"} if factor
is categorical, \code{c(value1, value2, \dots)}, or a function which
returns a vector, such as \code{seq(low,high,by=increment)}.  Only the
first factor may have the values omitted.  In this case the \code{Low
effect}, \code{Adjust to}, and \code{High effect} values will be used
from \code{datadist} if the variable is continuous.  For variables not
defined to \code{datadist}, you must specify non-missing constant
settings (or a vector of settings for the one displayed variable).  Note
that since \code{npsurv} objects do not use the variable list in 
\code{\dots}, you can specify any extra arguments to \code{labcurve} by
adding them at the end of the list of arguments.
}
\item{xlim}{
a vector of two numbers specifiying the x-axis range for follow-up time.
Default is \code{(0,maxtime)} where \code{maxtime} was the \code{pretty()}d version
of the maximum follow-up time
in any stratum, stored in \code{fit$maxtime}.  If \code{logt=TRUE},
default is \code{(1, log(maxtime))}.
}
\item{ylim}{
y-axis limits.  Default is \code{c(0,1)} for survival, and
\code{c(-5,1.5)} if \code{loglog=TRUE}. If \code{fun} or
\code{loglog=TRUE} are given and \code{ylim} is not, the limits will be
computed from the data.  For \code{what="hazard"}, default 
limits are computed from the first hazard function plotted.
}
\item{xlab}{
x-axis label.  Default is \code{units} attribute of failure time
variable given to \code{Surv}.
}
\item{ylab}{
y-axis label.  Default is \code{"Survival Probability"} or 
\code{"log(-log Survival Probability)"}. If \code{fun} is given, the default
is \code{""}.  For \code{what="hazard"}, the default is \code{"Hazard Function"}.
}
\item{time.inc}{
time increment for labeling the x-axis and printing numbers at risk. 
If not specified, the value
of \code{time.inc} stored with the model fit will be used.
}
\item{type}{
specifies type of estimates, \code{"tsiatis"} (the default) or \code{"kaplan-meier"}.
\code{"tsiatis"} here corresponds to the Breslow
estimator. This is ignored if survival estimates stored with \code{surv=TRUE} are
being used. For fits from \code{npsurv}, this argument
is also ignored, since it is specified as an argument to \code{npsurv}.
}
\item{conf.type}{
specifies the basis for confidence limits. This argument is
ignored for fits from \code{npsurv}.
}
\item{conf.int}{
Default is \code{FALSE}.  Specify e.g. \code{.95} to plot 0.95 confidence bands.
For fits from parametric survival models, or Cox models with
\code{x=TRUE} and \code{y=TRUE} specified to the fit, the exact
asymptotic formulas will be used to compute standard errors, and
confidence limits are based on \code{log(-log S(t))} if \code{loglog=TRUE}. 
If \code{x=TRUE} and \code{y=TRUE} were not specified to \code{cph} but
\code{surv=TRUE} was, the standard errors stored for the underlying
survival curve(s) will be used. These agree with the former if
predictions are requested at the mean value of X beta or if there are
only stratification factors in the model. This argument is ignored for
fits from \code{npsurv}, which must have previously specified
confidence interval specifications.  For \code{survdiffplot} if
\code{conf.int} is not specified, the level used in the call to
\code{npsurv} will be used.
}
\item{conf}{
\code{"bars"} for confidence bars at each \code{time.inc} time point. If the fit
was from \code{cph(\dots, surv=TRUE)}, the \code{time.inc} used will be
that stored with the fit. Use \code{conf="bands"} (the default) for
bands using standard errors at each failure time. For \code{npsurv}
objects only, \code{conf} may also be \code{"none"}, indicating that 
confidence interval information stored with the \code{npsurv} result
should be ignored.  For \code{npsurv} and \code{survdiffplot},
\code{conf} may be \code{"diffbands"} whereby a shaded region is drawn
for comparing two curves.  The polygon is centered at the midpoint of
the two survival estimates and the height of the polygon is 1/2 the
width of the approximate \code{conf.int} pointwise confidence region.
Survival curves not overlapping the shaded area are approximately
significantly different at the \code{1 - conf.int} level.
}
\item{what}{
defaults to \code{"survival"} to plot survival estimates.  Set to
\code{"hazard"} or an abbreviation to plot the hazard function (for
\code{psm} fits only). Confidence intervals are not available for
\code{what="hazard"}. 
}
\item{add}{
set to \code{TRUE} to add curves to an existing plot.
}
\item{label.curves}{
default is \code{TRUE} to use \code{labcurve} to label curves where they are farthest
apart.  Set \code{label.curves} to a \code{list} to specify options to
\code{labcurve}, e.g., \code{label.curves=list(method="arrow", cex=.8)}.
These option names may be abbreviated in the usual way arguments
are abbreviated.  Use for example \code{label.curves=list(keys=1:5)}
to draw symbols (as in \code{pch=1:5} - see \code{points})
on the curves and automatically position a legend
in the most empty part of the plot.  Set \code{label.curves=FALSE} to
suppress drawing curve labels.  The \code{col}, \code{lty}, \code{lwd}, and \code{type}
parameters are automatically passed to \code{labcurve}, although you
can override them here.  To distinguish curves by line types and
still have \code{labcurve} construct a legend, use for example
\code{label.curves=list(keys="lines")}.  The negative value for the
plotting symbol will suppress a plotting symbol from being drawn
either on the curves or in the legend.
}
\item{abbrev.label}{
set to \code{TRUE} to \code{abbreviate()} curve labels that are plotted
}
\item{levels.only}{
  set to \code{TRUE} to remove \code{variablename=} from the start of
  curve labels.}
\item{lty}{
vector of line types to use for different factor levels.  Default is
\code{c(1,3,4,5,6,7,\dots)}.
}
\item{lwd}{
vector of line widths to use for different factor levels.  Default is
current \code{par} setting for \code{lwd}.
}
\item{col}{
color for curve, default is \code{1}.  Specify a vector to assign different
colors to different curves.
}
\item{col.fill}{a vector of colors to used in filling confidence bands}
\item{adj.subtitle}{
set to \code{FALSE} to suppress plotting subtitle with levels of adjustment factors
not plotted. Defaults to \code{TRUE}.
This argument is ignored for \code{npsurv}.
}
\item{loglog}{
set to \code{TRUE} to plot \code{log(-log Survival)} instead of \code{Survival}
}
\item{fun}{
specifies any function to translate estimates and confidence limits
before plotting.
}
\item{logt}{
set to \code{TRUE} to plot \code{log(t)} instead of \code{t} on the x-axis
}
\item{n.risk}{
set to \code{TRUE} to add number of subjects at risk for each curve, using the
\code{surv.summary} created by \code{cph} or using the failure times used in
fitting the model if \code{y=TRUE} was specified to the fit or if the fit
was from \code{npsurv}.
The numbers are placed at the bottom
of the graph unless \code{y.n.risk} is given. 
If the fit is from \code{survest.psm}, \code{n.risk} does not apply.
}
\item{srt.n.risk}{
angle of rotation for leftmost number of subjects at risk (since this number
may run into the second or into the y-axis).  Default is \code{0}.
}
\item{adj.n.risk}{
justification for leftmost number at risk. Default is \code{1} for right 
justification.
Use \code{0} for left justification, \code{.5} for centered.
}
\item{sep.n.risk}{
multiple of upper y limit - lower y limit for separating lines of text
containing number of subjects at risk.  Default is \code{.056*(ylim[2]-ylim[1])}.
}
\item{y.n.risk}{
When \code{n.risk=TRUE}, the default is to place numbers of patients at risk above
the x-axis.  You can specify a y-coordinate for the bottom line of the
numbers using \code{y.n.risk}.
}
\item{cex.n.risk}{
character size for number of subjects at risk (when \code{n.risk} is \code{TRUE})
}
\item{dots}{
 set to \code{TRUE} to plot a grid of dots.  Will be plotted at every
 \code{time.inc} (see \code{cph}) and at survival increments of .1 (if
 \code{d>.4}), .05 (if \code{.2 < d <= .4}), or .025 (if \code{d <= .2}),
 where \code{d} is the range of survival displayed. 
}
\item{dotsize}{size of dots in inches}
\item{grid}{
 defaults to \code{NULL} (not drawing grid lines). Set to \code{TRUE} to
 plot \code{gray(.8)} grid lines, or specify any color.
}
\item{pr}{
set to \code{TRUE} to print survival curve coordinates used in the plots
}
\item{aehaz}{set to \code{TRUE} to add number of events and exponential
	distribution hazard rate estimates in curve labels}
\item{times}{a numeric vector of times at which to compute cumulative
	incidence probability estimates to add to curve labels}
\item{order}{
  an integer vector of length two specifying the order of groups when
  computing survival differences.  The default of \code{1:2} indicates
  that the second group is subtracted from the first.  Specify
  \code{order=2:1} to instead subtract the first from the second.  A
  subtitle indicates what was done.
}
}
\value{
list with components adjust (text string specifying adjustment levels)
and \code{curve.labels} (vector of text strings corresponding to levels of factor
used to distinguish curves). For \code{npsurv}, the returned value is the
vector of strata labels, or NULL if there are no strata.
}
\section{Side Effects}{
plots. If \code{par()$mar[4]<4}, issues \code{par(mar=)} to increment \code{mar[4]} by 2
if \code{n.risk=TRUE} and \code{add=FALSE}. The user may want to reset \code{par(mar)} in
this case to not leave such a wide right margin for plots. You usually
would issue \code{par(mar=c(5,4,4,2)+.1)}.
}
\details{
\code{survplot} will not work for Cox models with time-dependent covariables.
Use \code{survest} or \code{survfit} for that purpose.

There is a set a system option \code{\link[Hmisc]{mgp.axis.labels}} to allow x
and y-axes to have differing \code{mgp} graphical parameters (see \code{par}).
This is important when labels for y-axis tick marks are to be written
horizontally (\code{par(las=1)}), as a larger gap between the labels and
the tick marks are needed.  You can set the axis-specific 2nd
component of \code{mgp} using \code{mgp.axis.labels(c(xvalue,yvalue))}.
}
\seealso{
\code{\link{datadist}}, \code{\link{rms}}, \code{\link{cph}},
\code{\link{psm}}, \code{\link{survest}}, \code{\link{predictrms}},
\code{\link{plot.Predict}}, \code{\link{ggplot.Predict}},
\code{\link[Hmisc]{units}}, \code{\link[Hmisc]{errbar}},  
\code{\link{survfit}}, \code{\link[survival]{survreg.distributions}},
\code{\link[Hmisc]{labcurve}},
\code{\link[Hmisc]{mgp.axis}}, \code{\link{par}},
}
\examples{
# Simulate data from a population model in which the log hazard
# function is linear in age and there is no age x sex interaction
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('male','female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
dd <- datadist(age, sex)
options(datadist='dd')
S <- Surv(dt,e)

# When age is in the model by itself and we predict at the mean age,
# approximate confidence intervals are ok

f <- cph(S ~ age, surv=TRUE)
survplot(f, age=mean(age), conf.int=.95)
g <- cph(S ~ age, x=TRUE, y=TRUE)
survplot(g, age=mean(age), conf.int=.95, add=TRUE, col='red', conf='bars')

# Repeat for an age far from the mean; not ok
survplot(f, age=75, conf.int=.95)
survplot(g, age=75, conf.int=.95, add=TRUE, col='red', conf='bars')


#Plot stratified survival curves by sex, adj for quadratic age effect
# with age x sex interaction (2 d.f. interaction)

f <- cph(S ~ pol(age,2)*strat(sex), x=TRUE, y=TRUE)
#or f <- psm(S ~ pol(age,2)*sex)
Predict(f, sex, age=c(30,50,70))
survplot(f, sex, n.risk=TRUE, levels.only=TRUE)   #Adjust age to median
survplot(f, sex, logt=TRUE, loglog=TRUE)   #Check for Weibull-ness (linearity)
survplot(f, sex=c("male","female"), age=50)
                                        #Would have worked without datadist
                                        #or with an incomplete datadist
survplot(f, sex, label.curves=list(keys=c(2,0), point.inc=2))
                                        #Identify curves with symbols


survplot(f, sex, label.curves=list(keys=c('m','f')))
                                        #Identify curves with single letters


#Plots by quintiles of age, adjusting sex to male
options(digits=3)
survplot(f, age=quantile(age,(1:4)/5), sex="male")


#Plot survival Kaplan-Meier survival estimates for males
f <- npsurv(S ~ 1, subset=sex=="male")
survplot(f)


#Plot survival for both sexes
f <- npsurv(S ~ sex)
survplot(f)
#Check for log-normal and log-logistic fits
survplot(f, fun=qnorm, ylab="Inverse Normal Transform")
survplot(f, fun=function(y)log(y/(1-y)), ylab="Logit S(t)")

#Plot the difference between sexes
survdiffplot(f)

#Similar but show half-width of confidence intervals centered
#at average of two survival estimates
survplot(f, conf='diffbands')

options(datadist=NULL)
}
\keyword{survival}
\keyword{hplot}
\keyword{nonparametric}
\keyword{models}


