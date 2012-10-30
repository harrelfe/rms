\name{predab.resample}
\alias{predab.resample}
\title{Predictive Ability using Resampling}

\description{
\code{predab.resample} is a general-purpose
function that is used by functions for specific models.
It computes estimates of optimism of, and bias-corrected estimates of a vector
of indexes of predictive accuracy, for a model with a specified
design matrix, with or without fast backward step-down of predictors. If \code{bw=TRUE}, the design
matrix \code{x} must have been created by \code{ols}, \code{lrm}, or \code{cph}.
If \code{bw=TRUE}, \code{predab.resample} stores as the \code{kept}
attribute a logical matrix encoding which
factors were selected at each repetition.
}
\usage{
predab.resample(fit.orig, fit, measure, 
                method=c("boot","crossvalidation",".632","randomization"),
                bw=FALSE, B=50, pr=FALSE,
                rule="aic", type="residual", sls=.05, aics=0,
                tol=1e-12, force=NULL, estimates=TRUE,
                non.slopes.in.x=TRUE, kint=1,
                cluster, subset, group=NULL, debug=FALSE, \dots)
}
\arguments{
\item{fit.orig}{
object containing the original full-sample fit, with the \code{x=TRUE} and
\code{y=TRUE} options specified to the model fitting function.  This model
should be the FULL model including all candidate variables ever excluded
because of poor associations with the response.
}
\item{fit}{
a function to fit the model, either the original model fit, or a fit in a
sample.  fit has as arguments \code{x},\code{y}, \code{iter}, \code{penalty}, \code{penalty.matrix},
\code{xcol}, and other arguments passed to \code{predab.resample}. 
If you don't want \code{iter}
as an argument inside the definition of \code{fit}, add \dots to the end of its
argument list. \code{iter} is passed to \code{fit} to inform the function of the
sampling repetition number (0=original sample).  If \code{bw=TRUE}, \code{fit} should
allow for the possibility of selecting no predictors, i.e., it should fit an
intercept-only model if the model has intercept(s). \code{fit} must return
objects \code{coef} and \code{fail} (\code{fail=TRUE} if \code{fit} failed due to singularity or
non-convergence - these cases are excluded from summary statistics). \code{fit}
must add design attributes to the returned object if \code{bw=TRUE}.  
The \code{penalty.matrix} parameter is not used if \code{penalty=0}.  The \code{xcol}
vector is a vector of columns of \code{X} to be used in the current model fit.
For \code{ols} and \code{psm} it includes a \code{1} for the intercept position.
\code{xcol} is not defined if \code{iter=0} unless the initial fit had been from
a backward step-down.  \code{xcol} is used to select the correct rows and columns
of \code{penalty.matrix} for the current variables selected, for example.
}
\item{measure}{
a function to compute a vector of indexes of predictive accuracy for a given fit.
For \code{method=".632"} or \code{method="crossval"}, it will make the most sense for
measure to compute only indexes that are independent of sample size. The
measure function should take the following arguments or use \dots: \code{xbeta} 
(X beta for
current fit), \code{y}, \code{evalfit}, \code{fit}, \code{iter}, and \code{fit.orig}. \code{iter} is as in \code{fit}.
\code{evalfit} is set to \code{TRUE}
by \code{predab.resample} if the fit is being evaluated on the sample used to make the
fit, \code{FALSE} otherwise; \code{fit.orig} is the fit object returned by the original fit on the whole
sample. Using \code{evalfit} will sometimes save computations. For example, in
bootstrapping the area under an ROC curve for a logistic regression model,
\code{lrm} already computes the area if the fit is on the training sample. 
\code{fit.orig}
is used to pass computed configuration parameters from the original fit such as
quantiles of predicted probabilities that are used as cut points in other samples.
The vector created by measure should have \code{names()} associated with it.
}
\item{method}{
The default is \code{"boot"} for ordinary bootstrapping (Efron, 1983, Eq. 2.10).  
Use \code{".632"} for Efron's \code{.632} method (Efron, 1983, Section 6 and Eq. 6.10),
\code{"crossvalidation"}
for grouped cross--validation, \code{"randomization"} for the randomization method. May
be abbreviated down to any level, e.g. \code{"b"}, \code{"."}, \code{"cross"}, \code{"rand"}.
}
\item{bw}{
Set to \code{TRUE} to do fast backward step-down for each training sample. Default is \code{FALSE}.
}
\item{B}{
Number of repetitions, default=50. For \code{method="crossvalidation"}, this is also
the number of groups the original sample is split into.
}
\item{pr}{
\code{TRUE} to print results for each sample. Default is \code{FALSE}.
}
\item{rule}{
Stopping rule for fastbw, \code{"aic"} or \code{"p"}. Default is \code{"aic"} to use Akaike's
information criterion.
}
\item{type}{
Type of statistic to use in stopping rule for fastbw, \code{"residual"} (the default) or
\code{"individual"}.
}
\item{sls}{
Significance level for stopping in fastbw if \code{rule="p"}. Default is \code{.05}.
}
\item{aics}{
Stopping criteria for \code{rule="aic"}. Stops deleting factors when
chi-square - 2 times d.f. falls below \code{aics}. Default is \code{0}.
}
\item{tol}{
Tolerance for singularity checking.  Is passed to \code{fit} and \code{fastbw}.
}
\item{force}{see \code{\link{fastbw}}}
\item{estimates}{see \code{\link{print.fastbw}}}
\item{non.slopes.in.x}{set to \code{FALSE} if the design matrix \code{x}
does not have columns for intercepts and these columns are needed}
\item{kint}{
For multiple intercept models such as the ordinal logistic model, you may
specify which intercept to use as \code{kint}.  This affects the linear
predictor that is passed to \code{measure}.
}
\item{cluster}{
Vector containing cluster identifiers.  This can be specified only if
\code{method="boot"}.  If it is present, the bootstrap is done using sampling
with replacement from the clusters rather than from the original records.
If this vector is not the same length as the number of rows in the data
matrix used in the fit, an attempt will be made to use \code{naresid} on 
\code{fit.orig} to conform \code{cluster} to the data.  
See \code{bootcov} for more about this.
}
\item{subset}{
specify a vector of positive or negative integers or a logical vector when
you want to have the \code{measure} function compute measures of accuracy on
a subset of the data.  The whole dataset is still used for all model development.
For example, you may want to \code{validate} or \code{calibrate} a model by
assessing the predictions on females when the fit was based on males and
females.  When you use \code{cr.setup} to build extra observations for fitting the
continuation ratio ordinal logistic model, you can use \code{subset} to specify
which \code{cohort} or observations to use for deriving indexes of predictive
accuracy.  For example, specify \code{subset=cohort=="all"} to validate the
model for the first layer of the continuation ratio model (Prob(Y=0)).
}
\item{group}{
a grouping variable used to stratify the sample upon bootstrapping.
This allows one to handle k-sample problems, i.e., each bootstrap
sample will be forced to selected the same number of observations from
each level of group as the number appearing in the original dataset.
}
\item{debug}{set to \code{TRUE} to print subscripts of all training and
  test samples}
\item{\dots}{
The user may add other arguments here that are passed to \code{fit} and
\code{measure}.
}}
\value{
a matrix of class \code{"validate"} with rows corresponding
to indexes computed by \code{measure}, and the following columns:

\item{index.orig}{
indexes in original overall fit
}
\item{training}{
average indexes in training samples
}
\item{test}{
average indexes in test samples
}
\item{optimism}{
average \code{training-test} except for \code{method=".632"} - is .632 times
\code{(index.orig - test)}
}
\item{index.corrected}{
\code{index.orig-optimism}
}
\item{n}{
number of successful repetitions with the given index non-missing
}.
Also contains an attribute \code{keepinfo} if \code{measure} returned
such an attribute when run on the original fit.
}
\details{
For \code{method=".632"}, the program stops with an error if every observation
is not omitted at least once from a bootstrap sample.  Efron's ".632" method
was developed for measures that are formulated in terms on per-observation
contributions.  In general, error measures (e.g., ROC areas) cannot be
written in this way, so this function uses a heuristic extension to
Efron's formulation in which it is assumed that the average error measure
omitting the \code{i}th observation is the same as the average error measure
omitting any other observation.  Then weights are derived
for each bootstrap repetition and weighted averages over the \code{B} repetitions
can easily be computed.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Efron B, Tibshirani R (1997). Improvements on cross-validation: The .632+ bootstrap method.  JASA 92:548--560.
}
\seealso{
\code{\link{rms}}, \code{\link{validate}}, \code{\link{fastbw}},
\code{\link{lrm}}, \code{\link{ols}}, \code{\link{cph}},
\code{\link{bootcov}}
}
\examples{
# See the code for validate.ols for an example of the use of
# predab.resample
}
\keyword{models}
\concept{model validation}
\concept{bootstrap}
\concept{predictive accuracy}
