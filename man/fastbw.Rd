\name{fastbw}
\alias{fastbw}
\alias{print.fastbw}
\title{Fast Backward Variable Selection}
\description{
Performs a slightly inefficient but numerically stable version of fast
backward elimination on factors, using a method based on Lawless and Singhal
(1978).
This method uses the fitted complete model and computes approximate Wald
statistics by computing conditional (restricted) maximum likelihood estimates
assuming multivariate normality of estimates.
\code{fastbw} deletes factors, not columns of the design matrix. Factors requiring multiple d.f. will be retained or dropped as a group.
The function prints the deletion statistics for each variable in
turn, and prints approximate parameter estimates for the model after
deleting variables.  The approximation is better when the number of
factors deleted is not large.  For \code{ols}, the approximation is exact for
regression coefficients, and standard errors are only off by a factor
equal to the ratio of the mean squared error estimate for the reduced
model to the original mean squared error estimate for the full model.

If the fit was from \code{ols}, \code{fastbw} will compute the usual \eqn{R^2}
statistic for each model.
}
\usage{
fastbw(fit, rule="aic", type="residual", sls=.05, aics=0, eps=1e-9,
       k.aic=2, force=NULL)

\method{print}{fastbw}(x, digits=4, estimates=TRUE, \dots)
}
\arguments{
\item{fit}{
fit object with \code{Varcov(fit)} defined (e.g., from \code{ols}, \code{lrm}, \code{cph}, \code{psm}, \code{glmD})
}
\item{rule}{
Stopping rule. Defaults to \code{"aic"} for Akaike's information criterion. Use
\code{rule="p"} to use \eqn{P}-values
}
\item{type}{
Type of statistic on which to base the stopping rule. Default is
\code{"residual"} for
the pooled residual chi-square. Use \code{type="individual"} to use Wald
chi-square of individual factors.
}
\item{sls}{
Significance level for staying in a model if \code{rule="p"}.  Default is .05.
}
\item{aics}{
For \code{rule="aic"},
variables are deleted until the chi-square - \code{k.aic} times d.f. falls below \code{aics}.
Default \code{aics} is zero to use the ordinary AIC.  Set \code{aics} to say 10000
to see all variables deleted in order of descending importance.
}
\item{eps}{
Singularity criterion, default is \code{1E-9}.
}
\item{k.aic}{
multiplier to compute AIC, default is 2.  To use BIC, set \code{k.aic} equal
to \eqn{\log(n)}, where \eqn{n} is the effective sample size (number of events
for survival models).
}
\item{force}{a vector of integers specifying parameters forced to be in
  the model, not counting intercept(s)}
\item{x}{result of \code{fastbw}}
\item{digits}{number of significant digits to print}
\item{estimates}{set to \code{FALSE} to suppress printing table of
  approximate coefficients, SEs, etc., after variable deletions}
\item{\dots}{ignored}
}
\value{
a list with an attribute \code{kept} if \code{bw=TRUE}, and the
following components:

\item{result}{
matrix of statistics with rows in order of deletion.
}
\item{names.kept}{
names of factors kept in final model.
}
\item{factors.kept}{
the subscripts of factors kept in the final model
}
\item{factors.deleted}{
opposite of \code{factors.kept}.
}
\item{parms.kept}{
column numbers in design matrix corresponding to parameters kept in
the final model.
}
\item{parms.deleted}{
opposite of \code{parms.kept}.
}
\item{coefficients}{
vector of approximate coefficients of reduced model.
}
\item{var}{
approximate covariance matrix for reduced model.
}
\item{Coefficients}{
matrix of coefficients of all models.  Rows correspond to the
successive models examined and columns correspond to the coefficients
in the full model.  For variables not in a particular sub-model (row),
the coefficients are zero.
}}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Lawless, J. F. and Singhal, K. (1978): Efficient screening of nonnormal regression models.  Biometrics 34:318--327.
}
\seealso{
  \code{\link{rms}}, \code{\link{ols}}, \code{\link{lrm}},
  \code{\link{cph}}, \code{\link{psm}}, \code{\link{validate}},
  \code{\link[Hmisc]{solvet}}, \code{\link{rmsMisc}}
}
\examples{
\dontrun{
fastbw(fit, optional.arguments)     # print results
z <- fastbw(fit, optional.args)     # typically used in simulations
lm.fit(X[,z$parms.kept], Y)         # least squares fit of reduced model
}
}
\keyword{models}
\keyword{regression}
\keyword{htest}
\concept{stepwise}
\concept{variable selection}
