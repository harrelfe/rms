\name{val.surv}
\alias{val.surv}
\alias{plot.val.surv}
\alias{plot.val.survh}
\alias{print.val.survh}
\title{
Validate Predicted Probabilities Against Observed Survival Times
}
\description{
The \code{val.surv} function is useful for validating predicted survival
probabilities against right-censored failure times.  If \code{u} is
specified, the hazard regression function \code{hare} in the
\code{polspline} package is used to relate predicted survival
probability at time \code{u} to observed survival times (and censoring
indicators) to estimate the actual survival probability at time
\code{u} as a function of the estimated survival probability at that
time, \code{est.surv}.  If \code{est.surv} is not given, \code{fit} must
be specified and the \code{survest} function is used to obtain the
predicted values (using \code{newdata} if it is given, or using the
stored linear predictor values if not).  \code{hare} is given the sole
predictor \code{fun(est.surv)} where \code{fun} is given by the user or
is inferred from \code{fit}.  \code{fun} is the function of predicted
survival probabilities that one expects to create a linear relationship
with the linear predictors.

\code{hare} uses an adaptive procedure to find a linear spline of
\code{fun(est.surv)} in a model where the log hazard is a linear spline
in time \eqn{t}, and cross-products between the two splines are allowed so as to
not assume proportional hazards.  Thus \code{hare} assumes that the
covariate and time functions are smooth but not much else, if the number
of events in the dataset is large enough for obtaining a reliable
flexible fit.  There are special \code{print} and \code{plot} methods
when \code{u} is given.  In this case, \code{val.surv} returns an object
of class \code{"val.survh"}, otherwise it returns an object of class
\code{"val.surv"}.

If \code{u} is not specified, \code{val.surv} uses Cox-Snell (1968)
residuals on the cumulative 
probability scale to check on the calibration of a survival model
against right-censored failure time data.  If the predicted survival
probability at time \eqn{t} for a subject having predictors \eqn{X} is
\eqn{S(t|X)}, this method is based on the fact that the predicted
probability of failure before time \eqn{t}, \eqn{1 - S(t|X)}, when
evaluated at the subject's actual survival time \eqn{T}, has a uniform
(0,1) distribution.  The quantity \eqn{1 - S(T|X)} is right-censored
when \eqn{T} is.  By getting one minus the Kaplan-Meier estimate of the
distribution of \eqn{1 - S(T|X)} and plotting against the 45 degree line
we can check for calibration accuracy.  A more stringent assessment can
be obtained by stratifying this analysis by an important predictor
variable.  The theoretical uniform distribution is only an approximation
when the survival probabilities are estimates and not population values.

When \code{censor} is specified to \code{val.surv}, a different
validation is done that is more stringent but that only uses the
uncensored failure times.  This method is used for type I censoring when
the theoretical censoring times are known for subjects having uncensored
failure times.  Let \eqn{T}, \eqn{C}, and \eqn{F} denote respectively
the failure time, censoring time, and cumulative failure time
distribution (\eqn{1 - S}).  The expected value of \eqn{F(T | X)} is 0.5
when \eqn{T} represents the subject's actual failure time.  The expected
value for an uncensored time is the expected value of \eqn{F(T | T \leq
C, X) = 0.5 F(C | X)}.  A smooth plot of \eqn{F(T|X) - 0.5 F(C|X)} for
uncensored \eqn{T} should be a flat line through \eqn{y=0} if the model
is well calibrated.  A smooth plot of \eqn{2F(T|X)/F(C|X)} for
uncensored \eqn{T} should be a flat line through \eqn{y=1.0}. The smooth
plot is obtained by smoothing the (linear predictor, difference or
ratio) pairs.
}
\usage{
val.surv(fit, newdata, S, est.surv, censor,
         u, fun, lim, evaluate=100, pred, maxdim=5, ...)

\method{print}{val.survh}(x, ...)

\method{plot}{val.survh}(x, lim, xlab, ylab,
                         riskdist=TRUE, add=FALSE,
                         scat1d.opts=list(nhistSpike=200), ...)

\method{plot}{val.surv}(x, group, g.group=4,
     what=c('difference','ratio'),
     type=c('l','b','p'),
     xlab, ylab, xlim, ylim, datadensity=TRUE, \dots)
}
\arguments{
\item{fit}{a fit object created by \code{cph} or \code{psm}}
\item{newdata}{
a data frame for which \code{val.surv} should obtain predicted survival
probabilities.  If omitted, survival estimates are made for all of the
subjects used in \code{fit}.
}
\item{S}{an \code{\link[survival]{Surv}} object}
\item{est.surv}{
a vector of estimated survival probabilities corresponding to times in
the first column of \code{S}.
}
\item{censor}{
a vector of censoring times.  Only the censoring times for uncensored
observations are used.
}
\item{u}{a single numeric follow-up time}
\item{fun}{a function that transforms survival probabilities into the
     scale of the linear predictor.  If \code{fit} is given, and
     represents either a Cox, Weibull, or exponential fit, \code{fun} is
     automatically set to log(-log(p)).}
\item{lim}{a 2-vector specifying limits of predicted survival
     probabilities for obtaining estimated actual probabilities at time
     \code{u}.  Default for
     \code{val.surv} is the limits for predictions from \code{datadist},
     which for large \eqn{n} is the 10th smallest and 10th largest
     predicted survival probability.  For \code{plot.val.survh}, the
     default for \code{lim} is the range of the combination of predicted
     probabilities and calibrated actual probabilities.  \code{lim} is
     used for both axes of the calibration plot.}
\item{evaluate}{the number of evenly spaced points over the range of
     predicted probabilities.  This defines the points at which
     calibrated predictions are obtained for plotting.}
\item{pred}{a vector of points at which to evaluate predicted
     probabilities, overriding \code{lim}}
\item{maxdim}{see \code{\link[polspline]{hare}}}
\item{x}{result of \code{val.surv}}
\item{xlab}{x-axis label.  For \code{plot.survh}, defaults for
     \code{xlab} and \code{ylab} come from \code{u} and the units of
     measurement for the raw survival times.}
\item{ylab}{y-axis label}
\item{riskdist}{set to \code{FALSE} to not call \code{scat1d} to draw the
  distribution of predicted (uncalibrated) probabilities}
\item{add}{set to \code{TRUE} if adding to an existing plot}
\item{scat1d.opts}{a \code{list} of options to pass to \code{scat1d}.
     By default, the option \code{nhistSpike=200} is passed so that a spike
     histogram is used if the sample size exceeds 200.}
\item{\dots}{When \code{u} is given to \code{val.surv}, \dots represents
     optional arguments to \code{hare}.  It can represent arguments to
     pass to \code{plot} or \code{lines} for 
     \code{plot.val.survh}.  Otherwise, \dots contains optional
     arguments for \code{plsmo} or \code{plot}.   For
     \code{print.val.survh}, \dots is ignored.} 
\item{group}{
a grouping variable.  If numeric this variable is grouped into
\code{g.group} quantile groups (default is quartiles).  \code{group},
     \code{g.group}, \code{what}, and \code{type} apply when
     \code{u} is not given.}
\item{g.group}{
number of quantile groups to use when \code{group} is given and variable
is numeric.
}
\item{what}{
the quantity to plot when \code{censor} was in effect.  The default is to
show the difference between cumulative probabilities and their
expectation given the censoring time.  Set \code{what="ratio"} to show the
ratio instead.
}
\item{type}{
Set to the default (\code{"l"}) to plot the trend line only, \code{"b"}
to plot both individual subjects ratios and trend lines, or
\code{"p"} to plot only points. 
}
\item{xlim}{
}
\item{ylim}{
axis limits for \code{plot.val.surv} when the \code{censor} variable was used.
}
\item{datadensity}{
By default, \code{plot.val.surv} will show the data density on each curve
that is created as a result of \code{censor} being present.  Set
\code{datadensity=FALSE} to suppress these tick marks drawn by \code{scat1d}.
}
}
\value{a list of class \code{"val.surv"} or \code{"val.survh"}}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\references{
Cox DR, Snell EJ (1968):A general definition of residuals (with
discussion).  JRSSB 30:248--275.


Kooperberg C, Stone C, Truong Y (1995): Hazard regression.  JASA 90:78--94.


May M, Royston P, Egger M, Justice AC, Sterne JAC (2004):Development and
validation of a prognostic model for survival time data: application to
prognosis of HIV positive patients treated with antiretroviral therapy.
Stat in Med 23:2375--2398.


Stallard N (2009): Simple tests for th external validation of mortality
prediction scores.  Stat in Med 28:377--388.
}
\seealso{
\code{\link{validate}}, \code{\link{calibrate}}, \code{\link[polspline]{hare}},
\code{\link[Hmisc]{scat1d}}, \code{\link{cph}}, \code{\link{psm}},
\code{\link{groupkm}}
}
\examples{
# Generate failure times from an exponential distribution
set.seed(123)              # so can reproduce results
n <- 1000
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('Male','Female'), n, rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
units(t) <- 'Year'
label(t) <- 'Time to Event'
ev <- ifelse(t <= cens, 1, 0)
t <- pmin(t, cens)
S <- Surv(t, ev)

# First validate true model used to generate data

# If hare is available, make a smooth calibration plot for 1-year
# survival probability where we predict 1-year survival using the
# known true population survival probability
# In addition, use groupkm to show that grouping predictions into
# intervals and computing Kaplan-Meier estimates is not as accurate.

if('polspline' \%in\% row.names(installed.packages())) {
  s1 <- exp(-h*1)
  w <- val.surv(est.surv=s1, S=S, u=1,
                fun=function(p)log(-log(p)))
  plot(w, lim=c(.85,1), scat1d.opts=list(nhistSpike=200, side=1))
  groupkm(s1, S, m=100, u=1, pl=TRUE, add=TRUE)
}

# Now validate the true model using residuals

w <- val.surv(est.surv=exp(-h*t), S=S)
plot(w)
plot(w, group=sex)  # stratify by sex


# Now fit an exponential model and validate
# Note this is not really a validation as we're using the
# training data here
f <- psm(S ~ age + sex, dist='exponential', y=TRUE)
w <- val.surv(f)
plot(w, group=sex)


# We know the censoring time on every subject, so we can
# compare the predicted Pr[T <= observed T | T>c, X] to
# its expectation 0.5 Pr[T <= C | X] where C = censoring time
# We plot a ratio that should equal one
w <- val.surv(f, censor=cens)
plot(w)
plot(w, group=age, g=3)   # stratify by tertile of age
}
\keyword{models}
\keyword{regression}
\keyword{smooth}
\keyword{survival}
\concept{model validation}
\concept{predictive accuracy}
