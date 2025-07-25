\name{validate.lrm}
\alias{validate.lrm}
\alias{validate.orm}
\title{Resampling Validation of a Logistic or Ordinal Regression Model}
\description{
The \code{validate} function when used on an object created by
\code{lrm} or \code{orm} does resampling validation of a logistic
regression model, 
with or without backward step-down variable deletion.  It provides
bias-corrected Somers' \eqn{D_{xy}} rank correlation, R-squared index,
the intercept and slope of an overall logistic calibration equation, the
maximum absolute difference in predicted and calibrated probabilities
\eqn{E_{max}}, the discrimination index \eqn{D} (model L.R. \eqn{(\chi^2
- 1)/n}{(chi-square - 1)/n}), the unreliability index \eqn{U} =
difference in -2 log likelihood between un-calibrated \eqn{X\beta}{X
beta} and \eqn{X\beta}{X beta} with overall intercept and slope
calibrated to test sample / n, the overall quality index (logarithmic
probability score) \eqn{Q = D - U}, and the Brier or quadratic
probability score, \eqn{B} (the last 3 are not computed for ordinal
models), the \eqn{g}-index, and \code{gp}, the \eqn{g}-index on the
probability scale.  The corrected slope can be thought of as shrinkage
factor that takes into account overfitting.  For \code{orm} fits, a
subset of the above indexes is provided, Spearman's \eqn{\rho} is
substituted for \eqn{D_{xy}}, and a new index is reported: \code{pdm}, the mean
	absolute difference between 0.5 and the predicted probability that
	\eqn{Y\geq} the marginal median of \eqn{Y}.

See \code{\link{predab.resample}} for information about confidence limits.
}

\usage{
# fit <- lrm(formula=response ~ terms, x=TRUE, y=TRUE) or orm
\method{validate}{lrm}(fit, method="boot", B=40,
         bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0,
         force=NULL, estimates=TRUE,
         pr=FALSE,  kint, Dxy.method,
         emax.lim=c(0,1), \dots)
\method{validate}{orm}(fit, method="boot", B=40, bw=FALSE, rule="aic",
         type="residual",	sls=.05, aics=0, force=NULL, estimates=TRUE,
         pr=FALSE,  ...)
}
\arguments{
\item{fit}{
a fit derived by \code{lrm} or \code{orm}. The options \code{x=TRUE} and
\code{y=TRUE} must have been specified.
}
\item{method,B,bw,rule,type,sls,aics,force,estimates,pr}{see \code{\link{validate}} and \code{\link{predab.resample}}}
\item{kint}{
In the case of an ordinal model, specify which intercept to validate.
Default is the middle intercept.  For \code{validate.orm},
intercept-specific quantities are not validated so this does not matter.
}
\item{Dxy.method}{deprecated and ignored.  \code{lrm} through \code{lrm.fit} computes exact rank correlation coefficients as of version 6.9-0.}
\item{emax.lim}{
range of predicted probabilities over which to compute the maximum
error.  Default is entire range. 
}
\item{\dots}{
other arguments to pass to \code{lrm.fit} and to \code{predab.resample} (note especially
the \code{group}, \code{cluster}, and \code{subset} parameters)
}}
\value{
a matrix with rows corresponding to \eqn{D_{xy}},
\eqn{R^2}, \code{Intercept}, \code{Slope}, \eqn{E_{max}}, \eqn{D},
\eqn{U}, \eqn{Q}, \eqn{B}, \eqn{g}, \eqn{gp}, and
columns for the original index, resample estimates, indexes applied to
the whole or omitted sample using the model derived from the resample,
average optimism, corrected index, and number of successful re-samples.
For \code{validate.orm} not all columns are provided, Spearman's rho
is returned instead of \eqn{D_{xy}}, and \code{pdm} is reported.
}
\section{Side Effects}{
prints a summary, and optionally statistics for each re-fit
}
\details{
If the original fit was created using penalized maximum likelihood estimation,
the same \code{penalty.matrix} used with the original
fit are used during validation.

See \url{https://fharrell.com/post/bootcal/} for simulations of the accuracy of approximate bootstrap confidence intervals for overfitting-corrected Brier scores.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\references{
Miller ME, Hui SL, Tierney WM (1991): Validation techniques for
logistic regression models.  Stat in Med 10:1213--1226.


Harrell FE, Lee KL (1985):  A comparison of the
\emph{discrimination}
of discriminant analysis and logistic regression under multivariate
normality.  In Biostatistics: Statistics in Biomedical, Public Health,
and Environmental Sciences.  The Bernard G. Greenberg Volume, ed. PK
Sen. New York: North-Holland, p. 333--343.
}
\seealso{
\code{\link{predab.resample}}, \code{\link{fastbw}}, \code{\link{lrm}},
\code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{calibrate}}, 
\code{\link[Hmisc]{somers2}}, \code{\link{cr.setup}},
\code{\link{gIndex}}, \code{\link{orm}}
}
\examples{
n <- 1000    # define sample size
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))


# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


f <- lrm(y ~ sex*rcs(cholesterol)+pol(age,2)+blood.pressure, x=TRUE, y=TRUE)
#Validate full model fit
validate(f, B=10)              # normally B=300
validate(f, B=10, group=y)  
# two-sample validation: make resamples have same numbers of
# successes and failures as original sample


#Validate stepwise model with typical (not so good) stopping rule
validate(f, B=10, bw=TRUE, rule="p", sls=.1, type="individual")


\dontrun{
#Fit a continuation ratio model and validate it for the predicted
#probability that y=0
u <- cr.setup(y)
Y <- u$y
cohort <- u$cohort
attach(mydataframe[u$subs,])
f <- lrm(Y ~ cohort+rcs(age,4)*sex, penalty=list(interaction=2))
validate(f, cluster=u$subs, subset=cohort=='all') 
#see predab.resample for cluster and subset
}
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
\concept{model validation}
\concept{predictive accuracy}
\concept{bootstrap}
