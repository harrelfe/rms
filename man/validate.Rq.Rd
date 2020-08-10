\name{validate.Rq}
\alias{validate.Rq}
\title{Validation of a Quantile Regression Model}

\description{
The \code{validate} function when used on an object created by
\code{Rq} does resampling validation of a quantile regression
model, with or without backward step-down variable deletion.  Uses
resampling to estimate the optimism in various measures of predictive
accuracy which include mean absolute prediction error (MAD), Spearman
rho, the \eqn{g}-index, and the intercept and slope 
of an overall 
calibration \eqn{a + b\hat{y}}{a + b * (predicted y)}.  The "corrected"
slope can be thought of as shrinkage factor that takes into account
overfitting.  \code{validate.Rq} can also be used when a model for a
continuous response is going to be applied to a binary response. A
Somers' \eqn{D_{xy}} for this case is computed for each resample by
dichotomizing \code{y}. This can be used to obtain an ordinary receiver
operating characteristic curve area using the formula \eqn{0.5(D_{xy} +
1)}. See \code{predab.resample} for the list of
resampling methods.

The LaTeX \code{needspace} package must be in effect to use the
\code{latex} method. 
}
\usage{
# fit <- fitting.function(formula=response ~ terms, x=TRUE, y=TRUE)
\method{validate}{Rq}(fit, method="boot", B=40,
         bw=FALSE, rule="aic", type="residual", sls=0.05, aics=0, 
         force=NULL, estimates=TRUE, pr=FALSE, u=NULL, rel=">",
         tolerance=1e-7, \dots)
}
\arguments{
\item{fit}{
a fit derived by \code{Rq}. The options \code{x=TRUE} and \code{y=TRUE}
must have been specified.  See \code{validate} for a description of
arguments \code{method} - \code{pr}.
}
\item{method,B,bw,rule,type,sls,aics,force,estimates,pr}{see
  \code{\link{validate}} and \code{\link{predab.resample}} and
  \code{\link{fastbw}}}
\item{u}{
If specifed, \code{y} is also dichotomized at the cutoff \code{u} for
the purpose of getting a bias-corrected estimate of \eqn{D_{xy}}.
}
\item{rel}{
relationship for dichotomizing predicted \code{y}. Defaults to
\code{">"} to use \code{y>u}. \code{rel} can also be \code{"<"},
\code{">="}, and \code{"<="}. 
}
\item{tolerance}{
  ignored
}
\item{\dots}{
other arguments to pass to \code{predab.resample}, such as \code{group}, \code{cluster}, and \code{subset}
}}
\value{
matrix with rows corresponding to various indexes, and 
optionally \eqn{D_{xy}}, and
columns for the original index, resample estimates, 
indexes applied to whole or omitted sample using model derived from
resample, average optimism, corrected index, and number of successful resamples.
}
\section{Side Effects}{
prints a summary, and optionally statistics for each re-fit
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{Rq}}, \code{\link{predab.resample}}, \code{\link{fastbw}},
\code{\link{rms}}, \code{\link{rms.trans}},
\code{\link{gIndex}}
}
\examples{
set.seed(1)
x1 <- runif(200)
x2 <- sample(0:3, 200, TRUE)
x3 <- rnorm(200)
distance <- (x1 + x2/3 + rnorm(200))^2

f <- Rq(sqrt(distance) ~ rcs(x1,4) + scored(x2) + x3, x=TRUE, y=TRUE)

#Validate full model fit (from all observations) but for x1 < .75
validate(f, B=20, subset=x1 < .75)   # normally B=300

#Validate stepwise model with typical (not so good) stopping rule
validate(f, B=20, bw=TRUE, rule="p", sls=.1, type="individual")
}
\keyword{models}
\keyword{regression}
\concept{model validation}
\concept{bootstrap}
\concept{predictive accuracy}
