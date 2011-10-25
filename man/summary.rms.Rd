\name{summary.rms}
\alias{summary.rms}
\alias{print.summary.rms}
\alias{latex.summary.rms}
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
\code{print.summary.rms} prints the results, \code{latex.summary.rms} typesets
the results, and \code{plot.summary.rms}
plots shaded confidence bars to display the results graphically.
The longest confidence bar on each page is labeled with confidence levels
(unless this bar has been ignored due to \code{clip}).  By default, the
following confidence levels are all shown: .7, .8, .9, .95, and .99, using 
levels of gray scale (colors for Windows).
}
\usage{
\method{summary}{rms}(object, \dots, est.all=TRUE, antilog,
conf.int=.95, abbrev=FALSE, vnames=c("names","labels"))

\method{print}{summary.rms}(x, \dots)

\method{latex}{summary.rms}(object, title, table.env=TRUE, \dots)

\method{plot}{summary.rms}(x, at, log=FALSE,
    q=c(0.7, 0.8, 0.9, 0.95, 0.99), xlim, nbar, cex=1, nint=10,
    cex.c=.5, cex.t=1, clip=c(-1e30,1e30), main, \dots)
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
factors, a single value specifies the reference cell and the adjustment value. For
example, if \code{treat} has levels \code{"a", "b"} and \code{"c"} and \code{treat="b"}
is given to \code{summary}, treatment \code{a} will be compared to \code{b}
and \code{c} will be compared to \code{b}. Treatment \code{b} will be used when estimating
the effect of other factors. Category variables can have category labels
listed (in quotes), or an unquoted number that is a legal level, if all levels 
are numeric.  You need only use the first few
letters of each variable name - enough for unique identification.
For variables not defined with \code{datadist}, you must specify 3 values, none
of which are \code{NA}.

Also represents other arguments to pass to \code{latex}, is ignored for
\code{print}, or other optional arguments passed to \code{confbar}.  The
most important of these is \code{q}, the vector of confidence levels,
and \code{col}, which is a vector corresponding to \code{q} specifying
the colors for the regions of the bars.  \code{q} defaults to
\code{c(.7,.8,.9,.95,.99)} and \code{col} to \code{c(1,.8,.5,.2,.065)}
for UNIX, so that lower confidence levels (inner regions of bars)
corresponding with darker shades.  Specify for example \code{col=1:5} to
use actual colors.  For Windows, the default is \code{col=c(1,4,3,2,5)},
which by default represents black, blue, green, red, yellow.  If you
specify \code{q} you may have to specify \code{col}.
}
\item{est.all}{
Set to \code{FALSE} to only estimate effects of variables listed. Default is \code{TRUE}.
}
\item{antilog}{
Set to \code{FALSE} to suppress printing of anti-logged effects. Default is \code{TRUE}
if the model was fitted by \code{lrm} or \code{cph}.
Antilogged effects will be odds ratios for logistic models and hazard ratios
for proportional hazards models.
}
\item{conf.int}{
Defaults to \code{.95} for \code{95\%} confidence intervals of effects.
}
\item{abbrev}{
Set to \code{TRUE} to use the \code{abbreviate} function to shorten
factor levels for categorical variables in the model.
}
\item{vnames}{
  Set to \code{"labels"} to use variable labels to label effects.
  Default is \code{"names"} to use variable names.}
\item{x}{result of \code{summary}}
\item{title}{
\code{title} to pass to \code{latex}.  Default is name of fit object passed to
\code{summary} prefixed with \code{"summary"}.}
\item{table.env}{see \code{\link[Hmisc]{latex}}}
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
if \code{log=TRUE}).  If \code{at} is specified and \code{xlim} is omitted, \code{xlim} is
derived from the range of \code{at}.
}
\item{nbar}{
Sets up plot to leave room for \code{nbar} horizontal bars.  Default is the
number of non-interaction factors in the model.  Set \code{nbar} to a larger
value to keep too much surrounding space from appearing around horizontal
bars.  If \code{nbar} is smaller than the number of bars, the plot is divided
into multiple pages with up to \code{nbar} bars on each page.
}
\item{cex}{
\code{cex} parameter for factor labels.
}
\item{nint}{
Number of tick mark numbers for \code{pretty}.
}
\item{cex.c}{
\code{cex} parameter for \code{confbar}, for quantile labels.
}
\item{cex.t}{
\code{cex} parameter for main title.  Set to \code{0} to suppress the title.
}
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
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{datadist}}, \code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{rmsMisc}},
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
latex(s, file="")                # Just write LaTeX code to screen
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
