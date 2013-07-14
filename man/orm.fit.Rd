\name{orm.fit}
\alias{orm.fit}
\title{Ordinal Regression Model Fitter}
\description{
	Fits ordinal cumulative probability models for continuous or ordinal
	response variables, efficiently allowing for a large number of
	intercepts by capitalizing on the information matrix being sparse.
	Five different distribution functions are implemented, with the
	default being the logistic (yielding the proportional odds
	model).  Penalized estimation will be implemented in the future.
	Weights are not implemented.
  The optimization method is Newton-Raphson with step-halving.
}
\usage{
orm.fit(x=NULL, y, family='logistic',
        offset=0., initial,  maxit=12L, eps=.005, tol=1e-7, trace=FALSE,
        penalty.matrix=NULL)
}
\arguments{
\item{x}{
design matrix with no column for an intercept
}
\item{y}{
response vector, numeric, factor, or character.  The ordering of levels
is assumed from \code{factor(y)}.
}
\item{family}{the distribution family, corresponding to logistic (the
	default), Gaussian, Cauchy, Gumbel maximum (\eqn{exp(-exp(-x))};
	extreme value type I), and Gumbel minimum
	(\eqn{1-exp(-exp(x))}) distributions.  These are the cumulative
	distribution functions assumed for \eqn{Prob[Y \ge y | X]}.  The
	\code{family} argument can be an unquoted or a quoted string,
	e.g. \code{family=loglog} or \code{family="loglog"}.  To use
	a built-in family, the string must be one of the following
	corresponding to the previous list: \code{logistic, probit, loglog,
		cloglog, cauchit}.  The user can also provide her own customized
	family by setting \code{family} to a list with elements \code{cumprob,
	inverse, deriv, deriv2}; see the body of \code{orm.fit} for examples.
  An additional element, \code{name} must be given, which is a character
  string used to name the family for \code{print} and \code{latex}.}
\item{offset}{optional numeric vector containing an offset on the logit scale}
\item{initial}{vector of initial parameter estimates, beginning with the
	intercepts.  If \code{initial} is not specified, the function computes
  the overall score \eqn{\chi^2} test for the global null hypothesis of
	no regression.}
\item{maxit}{maximum no. iterations (default=\code{12}).}
\item{eps}{
difference in \eqn{-2 log} likelihood for declaring convergence.
Default is \code{.005}.  If the \eqn{-2 log} likelihood gets
worse by eps/10 while the maximum absolute first derivative of
\preformatted{-2 log} likelihood is below 1E-9, convergence is still
declared.  This handles the case where the initial estimates are MLEs,
to prevent endless step-halving.
}
\item{tol}{Singularity criterion. Default is 1e-7}
\item{trace}{
set to \code{TRUE} to print -2 log likelihood, step-halving
fraction, change in -2 log likelihood, and maximum absolute value of first
derivative at each iteration.
}
\item{penalty.matrix}{
	a self-contained ready-to-use penalty matrix - see\code{lrm} 
}
}

\value{
a list with the following components:

\item{call}{
calling expression
}
\item{freq}{
table of frequencies for \code{y} in order of increasing \code{y}
}
\item{yunique}{vector of sorted unique values of \code{y}}
\item{stats}{
vector with the following elements: number of observations used in the
fit, number of unique \code{y} values, median \code{y} from among the
observations used in the fit, maximum absolute value of first
derivative of log likelihood, model likelihood ratio chi-square, d.f.,
P-value, score chi-square and its P-value, Spearman's \eqn{\rho} rank
correlation between linear predictor and \code{y}, the
Nagelkerke \eqn{R^2} index, the \eqn{g}-index, \eqn{gr} (the
\eqn{g}-index on the ratio scale), and \eqn{pdm} (the mean absolute
difference between 0.5 and the estimated probability that \eqn{y\geq}
the marginal median).
When \code{penalty.matrix} is present, the \eqn{\chi^2}{chi-square},
d.f., and P-value are not corrected for the effective d.f.
}
\item{fail}{
set to \code{TRUE} if convergence failed (and \code{maxiter>1})
}
\item{coefficients}{
estimated parameters
}
\item{var}{
estimated variance-covariance matrix (inverse of information matrix).
Note that in the case of penalized estimation, \code{var} is not the
improved sandwich-type estimator (which \code{lrm} does compute).  The
only intercept parameter included in the stored object is the middle
intercept.
}
\item{family, trans}{see \code{\link{orm}}}
\item{deviance}{
-2 log likelihoods. 
When an offset variable is present, three
deviances are computed: for intercept(s) only, for
intercepts+offset, and for intercepts+offset+predictors.
When there is no offset variable, the vector contains deviances for
the intercept(s)-only model and the model with intercept(s) and predictors.
}
\item{non.slopes}{number of intercepts in model}
\item{interceptRef}{the index of the middle (median) intercept used in
	computing the linear predictor and \code{var}}
\item{linear.predictors}{the linear predictor using the first intercept}
\item{penalty.matrix}{see above}
\item{info.matrix}{see \code{\link{orm}}}
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
	\code{\link{orm}}, \code{\link{lrm}}, \code{\link{glm}},
	\code{\link{gIndex}}, \code{\link[SparseM]{solve}}
}
\examples{
#Fit an additive logistic model containing numeric predictors age, 
#blood.pressure, and sex, assumed to be already properly coded and 
#transformed
#
# fit <- orm.fit(cbind(age,blood.pressure,sex), death)
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
