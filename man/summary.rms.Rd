\name{summary.rms}
\alias{summary.rms}
\alias{print.summary.rms}
\alias{latex.summary.rms}
\alias{html.summary.rms}
\alias{plot.summary.rms}
\title{Summary of Effects in Model}

\description{
\code{summary.rms} forms a summary of the effects of each
factor.  When \code{summary} is used to estimate odds or hazard ratios for
continuous variables, it allows the levels of interacting factors to be
easily set, as well as allowing the user to choose the interval for the
effect. This method of estimating effects allows for nonlinearity in
the predictor.  Factors requiring multiple parameters are handled, as
\code{summary} obtains predicted values at the needed points and takes
differences.  By default, inter-quartile range effects (odds ratios,
hazards ratios, etc.) are printed for continuous factors, and all
comparisons with the reference level are made for categorical factors.
\code{print.summary.rms} prints the results, \code{latex.summary.rms}
and \code{html.summary.rms} typeset the results, and \code{plot.summary.rms}
plots shaded confidence bars to display the results graphically.
The longest confidence bar on each page is labeled with confidence levels
(unless this bar has been ignored due to \code{clip}).  By default, the
following confidence levels are all shown: .9, .95, and .99, using 
blue of different transparencies.  The \code{plot} method currently
ignores bootstrap and Bayesian highest posterior density intervals but approximates
intervals based on standard errors.  The \code{html} method is for use
with R Markdown using html.

The \code{print} method will call the \code{latex} or \code{html} method
if \code{options(prType=)} is set to \code{"latex"} or \code{"html"}.
For \code{"latex"} printing through \code{print()}, the LaTeX table
environment is turned off.

If \code{usebootcoef=TRUE} and the fit was run through \code{bootcov},
the confidence intervals are bootstrap nonparametric percentile
confidence intervals, basic bootstrap, or BCa intervals, obtained on contrasts
evaluated on all bootstrap samples.

If \code{options(grType='plotly')} is in effect and the \code{plotly}
package is installed, \code{plot} is used instead of base graphics to
draw the point estimates and confidence limits when the \code{plot}
method for \code{summary} is called.  Colors and other graphical
arguments to \code{plot.summary} are ignored in this case.  Various
special effects are implemented such as only drawing 0.95 confidence
limits by default but including a legend that allows the other CLs to be
activated.  Hovering over point estimates shows adjustment values if
there are any.  \code{nbar} is not implemented for \code{plotly}.
}
\usage{
\method{summary}{rms}(object, \dots, ycut=NULL, est.all=TRUE, antilog,
conf.int=.95, abbrev=FALSE, vnames=c("names","labels"),
conf.type=c('individual','simultaneous'),
usebootcoef=TRUE, boot.type=c("percentile","bca","basic"),
posterior.summary=c('mean', 'median', 'mode'), verbose=FALSE)

\method{print}{summary.rms}(x, \dots, table.env=FALSE)

\method{latex}{summary.rms}(object, title, table.env=TRUE, \dots)

\method{html}{summary.rms}(object, digits=4, dec=NULL, \dots)

\method{plot}{summary.rms}(x, at, log=FALSE,
    q=c(0.9, 0.95, 0.99), xlim, nbar, cex=1, nint=10,
    cex.main=1, clip=c(-1e30,1e30), main,
    col=rgb(red=.1,green=.1,blue=.8,alpha=c(.1,.4,.7)),
    col.points=rgb(red=.1,green=.1,blue=.8,alpha=1), pch=17,
    lwd=if(length(q) == 1) 3 else 2 : (length(q) + 1), digits=4, \dots)
}
\arguments{
\item{object}{
a \code{rms} fit object.  Either \code{options(datadist)} should have
been set before the fit, or \code{datadist()} and
\code{options(datadist)} run before \code{summary}.  For \code{latex} is
the result of \code{summary}.
}
\item{\dots}{
For \code{summary}, omit list of variables to estimate effects for all
predictors. Use a list 
of variables of the form \code{age, sex} to estimate using default
ranges. Specify \code{age=50} for example to adjust age to 50 when testing
other factors (this will only matter for factors that interact with age).
Specify e.g. \code{age=c(40,60)} to estimate the effect of increasing age from
40 to 60. Specify \code{age=c(40,50,60)} to let age range from 40 to 60 and
be adjusted to 50 when testing other interacting factors. For category
factors, a single value specifies the reference cell and the adjustment
    value. For example, if \code{treat} has levels \code{"a", "b"} and
    \code{"c"} and \code{treat="b"} is given to \code{summary},
    treatment \code{a} will be compared to \code{b} and \code{c} will be
    compared to \code{b}. Treatment \code{b} will be used when
    estimating the effect of other factors. Category variables can have
    category labels listed (in quotes), or an unquoted number that is a
    legal level, if all levels  are numeric.  You need only use the
    first few letters of each variable name - enough for unique
    identification. For variables not defined with \code{datadist}, you
    must specify 3 values, none of which are \code{NA}.

Also represents other arguments to pass to \code{latex}, is ignored for
\code{print} and \code{plot}.
}
\item{ycut}{must be specified if the fit is a partial proportional odds
	model.  Specifies the single value of the response variable used to
	estimate ycut-specific regression effects, e.g., odds ratios}
\item{est.all}{
Set to \code{FALSE} to only estimate effects of variables listed. Default is \code{TRUE}.
}
\item{antilog}{
Set to \code{FALSE} to suppress printing of anti-logged effects. Default
is \code{TRUE} if the model was fitted by \code{lrm} or \code{cph}.
Antilogged effects will be odds ratios for logistic models and hazard ratios
for proportional hazards models.
}
\item{conf.int}{
Defaults to \code{.95} for \code{95\%} confidence intervals of effects.}
\item{abbrev}{
Set to \code{TRUE} to use the \code{abbreviate} function to shorten
factor levels for categorical variables in the model.}
\item{vnames}{
  Set to \code{"labels"} to use variable labels to label effects.
  Default is \code{"names"} to use variable names.}
\item{conf.type}{
  The default type of confidence interval computed for a given
  individual (1 d.f.) contrast is a pointwise confidence interval.  Set
  \code{conf.type="simultaneous"} to use the \code{multcomp} package's
  \code{glht} and \code{confint} functions to compute confidence
  intervals with simultaneous (family-wise) coverage, thus adjusting for
  multiple comparisons.  Contrasts are simultaneous only over groups of
  intervals computed together.  
}
\item{usebootcoef}{
  If \code{fit} was the result of \code{bootcov} but you want to use the
  bootstrap covariance matrix instead of the nonparametric percentile,
  basic, or BCa methods for confidence intervals (which uses all the bootstrap
  coefficients), specify \code{usebootcoef=FALSE}.}
\item{boot.type}{set to \code{'bca'} to compute BCa confidence
	limits or to \code{'basic'} to use the basic bootstrap.  The default
	is to compute percentile intervals.}
\item{posterior.summary}{set to \code{'mode'} or \code{'median'} to use the posterior
	mean/median instead of the mean for point estimates of contrasts}
\item{verbose}{set to \code{TRUE} when \code{conf.type='simultaneous'}
  to get output describing scope of simultaneous adjustments}
\item{x}{result of \code{summary}}
\item{title}{
\code{title} to pass to \code{latex}.  Default is name of fit object passed to
\code{summary} prefixed with \code{"summary"}.}
\item{table.env}{see \code{\link[Hmisc]{latex}}}
\item{digits,dec}{for \code{html.summary.rms}; \code{digits} is the
	number of significant digits for printing for effects, standard
	errors, and confidence limits.  It is ignored if \code{dec} is
	given. The statistics are rounded to \code{dec} digits to the right of
	the decimal point of \code{dec} is given.  \code{digits} is also the
	number of significant digits to format numeric hover text and labels
	for \code{plotly}.}
\item{at}{
vector of coordinates at which to put tick mark labels on the main axis.  If
\code{log=TRUE}, \code{at} should be in anti-log units.
}
\item{log}{
Set to \code{TRUE} to plot on \eqn{X\beta}{X beta} scale but labeled with
anti-logs. 
}
\item{q}{scalar or vector of confidence coefficients to depict}
\item{xlim}{
X-axis limits for \code{plot} in units of the linear predictors (log scale
if \code{log=TRUE}).  If \code{at} is specified and \code{xlim} is
omitted, \code{xlim} is derived from the range of \code{at}.
}
\item{nbar}{
Sets up plot to leave room for \code{nbar} horizontal bars.  Default is the
number of non-interaction factors in the model.  Set \code{nbar} to a larger
value to keep too much surrounding space from appearing around horizontal
bars.  If \code{nbar} is smaller than the number of bars, the plot is divided
into multiple pages with up to \code{nbar} bars on each page.
}
\item{cex}{\code{cex} parameter for factor labels.}
\item{nint}{Number of tick mark numbers for \code{pretty}.}
\item{cex.main}{\code{cex} parameter for main title.  Set to \code{0} to
  suppress the title.}
\item{clip}{
confidence limits outside the interval \code{c(clip[1], clip[2])} will be
ignored, and \code{clip} also be respected when computing \code{xlim}
when \code{xlim} is not specified.  \code{clip} should be in the units of
\code{fun(x)}.  If \code{log=TRUE}, \code{clip} should be in \eqn{X\beta}{X
  beta} units. 
}
\item{main}{
main title.  Default is inferred from the model and value of \code{log},
e.g., \code{"log Odds Ratio"}.
}
\item{col}{vector of colors, one per value of \code{q}}
\item{col.points}{color for points estimates}
\item{pch}{symbol for point estimates.  Default is solid triangle.}
\item{lwd}{line width for confidence intervals, corresponding to
	\code{q}}
}
\value{
For \code{summary.rms}, a matrix of class \code{summary.rms} 
with rows corresponding to factors in
the model and columns containing the low and high values for the effects,
the range for the effects, the effect point estimates (difference in
predicted values for high and low factor values), the standard error
of this effect estimate, and the lower and upper confidence limits.
If \code{fit$scale.pred} has a second level, two rows appear for each factor,
the second corresponding to anti--logged effects. Non--categorical factors
are stored first, and effects for any categorical factors are stored at
the end of the returned matrix.  \code{scale.pred} and \code{adjust}.  \code{adjust}
is a character string containing levels of adjustment variables, if
there are any interactions.  Otherwise it is "".
\code{latex.summary.rms} returns an object of class \code{c("latex","file")}.
It requires the \code{latex} function in Hmisc.
}
\author{
Frank Harrell\cr
Hui Nian\cr  
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{datadist}}, \code{\link{rms}}, \code{\link{rms.trans}},
\code{\link{rmsMisc}}, 
\code{\link[Hmisc]{Misc}}, \code{\link{pretty}}, \code{\link{contrast.rms}}
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


fit <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)))


s <- summary(fit)                # Estimate effects using default ranges
                                 # Gets odds ratio for age=3rd quartile
                                 # compared to 1st quartile
\dontrun{
latex(s)                         # Use LaTeX to print nice version
latex(s, file="")                # Just write LaTeX code to console
html(s)                          # html/LaTeX to console for knitr
# Or:
options(prType='latex')
summary(fit)                     # prints with LaTeX, table.env=FALSE
options(prType='html')
summary(fit)                     # prints with html
}
summary(fit, sex='male', age=60) # Specify ref. cell and adjustment val
summary(fit, age=c(50,70))       # Estimate effect of increasing age from
                                 # 50 to 70
s <- summary(fit, age=c(50,60,70)) 
                                 # Increase age from 50 to 70, adjust to
                                 # 60 when estimating effects of other factors
#Could have omitted datadist if specified 3 values for all non-categorical
#variables (1 value for categorical ones - adjustment level)
plot(s, log=TRUE, at=c(.1,.5,1,1.5,2,4,8))


options(datadist=NULL)
}
\keyword{models}
\keyword{regression}
\keyword{htest}
\keyword{survival}
\keyword{hplot}
\keyword{interface}
\concept{logistic regression model}
