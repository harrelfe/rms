\name{orm.fit}
\alias{orm.fit}
\title{Ordinal Regression Model Fitter}
\description{
	Fits ordinal cumulative probability models for continuous or ordinal
	response variables, efficiently allowing for a large number of
	intercepts by capitalizing on the information matrix being sparse.
	Five different distribution functions are implemented, with the
	default being the logistic (yielding the proportional odds
	model).  Penalized estimation and weights are also implemented, as in `[lrm.fit()]`.
  The optimization method is Newton-Raphson with step-halving, or the Levenberg-Marquart method.
  The latter has been shown to converge better when there are large offsets.
	Execution time is is fast even for hundreds of thousands of intercepts.  The limiting factor
  is the number of intercepts times the number of columns of \code{x}.
}
\usage{
orm.fit(x=NULL, y, family=c("logistic","probit","loglog","cloglog","cauchit"),
        offset, initial, opt_method=c('NR', 'LM'),
        maxit=30L, eps=5e-4, gradtol=0.001, abstol=1e10, 
        minstepsize=0.01, tol=.Machine$double.eps, trace=FALSE,
        penalty.matrix=NULL, weights=NULL, normwt=FALSE, scale=FALSE, mscore=FALSE,
        inclpen=TRUE, y.precision = 7, compstats=TRUE, onlydata=FALSE, ...)
}
\arguments{
\item{x}{
design matrix with no column for an intercept
}
\item{y}{
response vector, numeric, factor, or character.  The ordering of levels
is assumed from \code{factor(y)}.
}
\item{family}{a character value specifying the distribution family, corresponding to logistic (the
	default), Gaussian, Cauchy, Gumbel maximum (\eqn{exp(-exp(-x))};
	extreme value type I), and Gumbel minimum
	(\eqn{1-exp(-exp(x))}) distributions.  These are the cumulative
	distribution functions assumed for \eqn{Prob[Y \ge y | X]}.  The
	\code{family} argument can be an unquoted or a quoted string,
	e.g. \code{family=loglog} or \code{family="loglog"}.  To use
	a built-in family, the string must be one of the following
	corresponding to the previous list: \code{logistic, probit, loglog,
		cloglog, cauchit}.}
\item{offset}{optional numeric vector containing an offset on the logit scale}
\item{initial}{vector of initial parameter estimates, beginning with the
	intercepts.  If \code{initial} is not specified, the function computes
  the overall score \eqn{\chi^2} test for the global null hypothesis of
	no regression.  \code{initial} is padded to the right with zeros for the regression
  coefficients, if needed.  When censoring is present, \code{initial} can also be
  a list with elements \code{time} and \code{surv} from the \code{npsurv} attribute of the
  \code{y} element of a previous fit.  This is useful when bootstrapping, for example.}
\item{opt_method}{set to \code{"LM"} to use Levenberg-Marquardt instead of the default Newton-Raphson}
\item{maxit}{maximum no. iterations (default=\code{30}).}
\item{eps}{
difference in \eqn{-2 log} likelihood for declaring convergence.
Default is \code{.0005}.  This handles the case where the initial estimates are MLEs,
to prevent endless step-halving.
}
\item{gradtol}{maximum absolute gradient before convergence can be declared. \code{gradtol} is automatically scaled by n / 1000 since the gradient is proportional to the sample size.}
\item{abstol}{maximum absolute change in parameter estimates from one iteration to the next before convergence can be declared; by default has no effect}
\item{minstepsize}{used to specify when to abandon step-halving}
\item{tol}{Singularity criterion. Default is typically 2e-16}
\item{trace}{
set to \code{TRUE} to print -2 log likelihood, step-halving
fraction, change in -2 log likelihood, maximum absolute value of first
derivative, and max absolute change in parameter estimates at each iteration.
}
\item{penalty.matrix}{
	a self-contained ready-to-use penalty matrix - see\code{lrm}
}
\item{weights}{a vector (same length as \code{y}) of possibly fractional case weights}
\item{normwt}{set to \code{TRUE} to scale \code{weights} so they sum to \eqn{n}, the length of \code{y}; useful for sample surveys as opposed to the default of frequency weighting}
\item{mscore}{set to \code{TRUE} to compute the sparse score matrix and store its elements as a list \code{mscore}}
\item{scale}{set to \code{TRUE} to subtract column means and divide by
	column standard deviations of \code{x}
  before fitting, and to back-solve for the un-normalized covariance
  matrix and regression coefficients.  This can sometimes make the model
	converge for very large
  sample sizes where for example spline or polynomial component
  variables create scaling problems leading to loss of precision when
  accumulating sums of squares and crossproducts.}
\item{inclpen}{set to \code{FALSE} to not include the penalty matrix in the Hessian when the Hessian is being computed on transformed \code{x}, vs. adding the penalty after back-transforming.  This should not matter.}
\item{y.precision}{When \sQuote{y} is numeric, values may need to be rounded
to avoid unpredictable behavior with \code{unique()} with floating-point
numbers. Default is to 7 decimal places.}
\item{compstats}{set to \code{FALSE} to prevent the calculation of the vector of model statistics}
\item{onlydata}{set to \code{TRUE} to return the data used in model fitting as a list, without fitting the model}
\item{\dots}{ignored}
}
\value{
a list with the following components, not counting all the components produced by `orm.fit`:

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
correlation between linear predictor and \code{y} (if there is no censoring),
Somers' \eqn{Dxy} rank correlation (if there is no censoring or only
right censoring),) the
Nagelkerke \eqn{R^2} index, other \eqn{R^2} measures, the \eqn{g}-index, \eqn{gr} (the
\eqn{g}-index on the ratio scale), and \eqn{pdm} (the mean absolute
difference between 0.5 and the estimated probability that \eqn{y\geq}
the marginal median).
When \code{penalty.matrix} is present, the \eqn{\chi^2}{chi-square},
d.f., and P-value are not corrected for the effective d.f.
}
\item{fail}{
set to \code{TRUE} if convergence failed (and \code{maxit>1})
}
\item{coefficients}{
estimated parameters
}
\item{family, famfunctions}{see \code{\link{orm}}}
\item{deviance}{
-2 log likelihoods.
When an offset variable is present, three
deviances are computed: for intercept(s) only, for
intercepts+offset, and for intercepts+offset+predictors.
When there is no offset variable, the vector contains deviances for
the intercept(s)-only model and the model with intercept(s) and predictors.
}
\item{lpe}{vector of per-observation likelihood probability elements.  An observation's
contribution to the log likelihood is the log of \code{lpe}.}
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
fh@fharrell.com
}
\seealso{
	\code{\link{orm}}, \code{\link{lrm}}, \code{\link{glm}},
	\code{\link{gIndex}}, \code{\link[SparseM:SparseM.solve]{solve}},
  \code{\link{recode2integer}}
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
