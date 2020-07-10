\name{orm}
\alias{orm}
\alias{print.orm}
\alias{Quantile.orm}
\title{Ordinal Regression Model}
\description{
	Fits ordinal cumulative probability models for continuous or ordinal
	response variables, efficiently allowing for a large number of
	intercepts by capitalizing on the information matrix being sparse.
	Five different distribution functions are implemented, with the
	default being the logistic (i.e., the proportional odds
	model).  The ordinal cumulative probability models are stated in terms
	of exceedance probabilities (\eqn{Prob[Y \ge y | X]}) so that as with
	OLS larger predicted values are associated with larger \code{Y}.  This is
	important to note for the asymmetric distributions given by the
	log-log and complementary log-log families, for which negating the
	linear predictor does not result in \eqn{Prob[Y < y | X]}.  The
	\code{family} argument is defined in \code{\link{orm.fit}}.  The model
	assumes that the inverse of the assumed cumulative distribution
	function, when applied to one minus the true cumulative distribution function
	and plotted on the \eqn{y}-axis (with the original \eqn{y} on the
	\eqn{x}-axis) yields parallel curves (though not necessarily linear).
	This can be checked by plotting the inverse cumulative probability
	function of one minus the empirical distribution function, stratified
	by \code{X}, and assessing parallelism.  Note that parametric
	regression models make the much stronger assumption of linearity of
	such inverse functions.

	For the \code{print} method, format of output is controlled by the
	user previously running \code{options(prType="lang")} where
	\code{lang} is \code{"plain"} (the default), \code{"latex"}, or
	\code{"html"}. 

	\code{Quantile.orm} creates an R function that computes an estimate of
	a given quantile for a given value of the linear predictor (using the
	first intercept).  It uses approximate linear interpolation on the
	original response variable scale (the interpolation is exactly linear
	of there is only one unique value of the linear predictor given as the
	second argument).
}
\usage{
orm(formula, data=environment(formula),
    subset, na.action=na.delete, method="orm.fit",
    model=FALSE, x=FALSE, y=FALSE, linear.predictors=TRUE, se.fit=FALSE, 
    penalty=0, penalty.matrix, tol=1e-7, eps=0.005, 
    var.penalty=c('simple','sandwich'), scale=FALSE, \dots)

\method{print}{orm}(x, digits=4, coefs=TRUE,
    intercepts=x$non.slopes < 10, title, \dots)

\method{Quantile}{orm}(object, codes=FALSE, \dots)
}
\arguments{
\item{formula}{
a formula object. An \code{offset} term can be included. The offset causes
fitting of a model such as \eqn{logit(Y=1) = X\beta + W}, where \eqn{W} is the
offset variable having no estimated coefficient.
The response variable can be any data type; \code{orm} converts it
in alphabetic or numeric order to a factor variable and
recodes it 1,2,\dots internally. 
}
\item{data}{
data frame to use. Default is the current frame.
}
\item{subset}{
logical expression or vector of subscripts defining a subset of
observations to analyze
}
\item{na.action}{
function to handle \code{NA}s in the data. Default is \code{na.delete}, which
deletes any observation having response or predictor missing, while
preserving the attributes of the predictors and maintaining frequencies
of deletions due to each variable in the model.  
This is usually specified using \code{options(na.action="na.delete")}.
}
\item{method}{
name of fitting function. Only allowable choice at present is \code{orm.fit}.
}
\item{model}{
causes the model frame to be returned in the fit object
}
\item{x}{
causes the expanded design matrix (with missings excluded)
to be returned under the name \code{x}.  For \code{print}, an object
created by \code{orm}.
}
\item{y}{
causes the response variable (with missings excluded) to be returned
under the name \code{y}.
}
\item{linear.predictors}{
causes the predicted X beta (with missings excluded) to be returned
under the name \code{linear.predictors}.  The first intercept is used.
}
\item{se.fit}{
causes the standard errors of the fitted values (on the linear predictor
scale) to be returned under the name \code{se.fit}.  The middle
intercept is used.
}
\item{penalty}{see \code{\link{lrm}}}
\item{penalty.matrix}{see \code{\link{lrm}}}
\item{tol}{singularity criterion (see \code{orm.fit})}
\item{eps}{difference in \eqn{-2 log} likelihood for declaring convergence}
\item{var.penalty}{see \code{\link{lrm}}}
\item{scale}{set to \code{TRUE} to subtract column means and divide by
	column standard deviations of the design matrix
  before fitting, and to back-solve for the un-normalized covariance
  matrix and regression coefficients.  This can sometimes make the model
	converge for very large 
  sample sizes where for example spline or polynomial component
  variables create scaling problems leading to loss of precision when
  accumulating sums of squares and crossproducts.}
\item{\dots}{arguments that are passed to \code{orm.fit}, or from
    \code{print}, to \code{\link{prModFit}}.  Ignored for
		\code{Quantile}.  One of the most important arguments is \code{family}.}
\item{digits}{number of significant digits to use}
\item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{intercepts}{By default, intercepts are only printed if there are
    fewer than 10 of them.  Otherwise this is controlled by specifying
    \code{intercepts=FALSE} or \code{TRUE}.}
\item{title}{a character string title to be passed to \code{prModFit}.
  Default is constructed from the name of the distribution family.}
\item{object}{an object created by \code{orm}}
\item{codes}{if \code{TRUE}, uses the integer codes \eqn{1,2,\ldots,k}
	for the \eqn{k}-level response in computing the predicted quantile}
}
\value{
The returned fit object of \code{orm} contains the following components
in addition to the ones mentioned under the optional arguments.

\item{call}{calling expression}
\item{freq}{
table of frequencies for \code{Y} in order of increasing \code{Y}}
\item{stats}{
vector with the following elements: number of observations used in the
fit, number of unique \code{Y} values, median \code{Y} from among the
observations used int he fit, maximum absolute value of first
derivative of log likelihood, model likelihood ratio
\eqn{\chi^2}{chi-square}, d.f.,  \eqn{P}-value, score \eqn{\chi^2}
statistic (if no initial values given), \eqn{P}-value, Spearman's
\eqn{\rho} rank correlation between the linear predictor and \code{Y},
the Nagelkerke \eqn{R^2} index, the \eqn{g}-index, \eqn{gr} (the
\eqn{g}-index on the odds ratio scale), and \eqn{pdm} (the mean absolute
difference between 0.5 and the predicted probability that \eqn{Y\geq}
the marginal median).
In the case of penalized estimation, the \code{"Model L.R."} is computed
without the penalty factor, and \code{"d.f."} is the effective d.f. from
Gray's (1992) Equation 2.9.  The \eqn{P}-value uses this corrected model
L.R. \eqn{\chi^2}{chi-square} and corrected d.f. 
The score chi-square statistic uses first derivatives which contain
penalty components.
}
\item{fail}{
set to \code{TRUE} if convergence failed (and \code{maxiter>1}) or if a
singular information matrix is encountered
}
\item{coefficients}{estimated parameters}
\item{var}{
estimated variance-covariance matrix (inverse of information matrix)
for the middle intercept and regression coefficients.  See
\code{\link{lrm}} for details if penalization is used.
}
\item{effective.df.diagonal}{see \code{\link{lrm}}}
\item{family}{the character string for \code{family}.  If \code{family}
	was a user-customized list, it must have had an element named
	\code{name}, which is taken as the return value for \code{family} here.}
\item{trans}{a list of functions for the choice of \code{family}, with
	elements \code{cumprob} (the cumulative probability distribution
	function), \code{inverse} (inverse of \code{cumprob}), \code{deriv}
	(first derivative of \code{cumprob}), and \code{deriv2} (second
	derivative of \code{cumprob})}
\item{deviance}{
-2 log likelihoods (counting penalty components)
When an offset variable is present, three
deviances are computed: for intercept(s) only, for
intercepts+offset, and for intercepts+offset+predictors.
When there is no offset variable, the vector contains deviances for
the intercept(s)-only model and the model with intercept(s) and predictors.
}
\item{non.slopes}{number of intercepts in model}
\item{interceptRef}{the index of the middle (median) intercept used in
	computing the linear predictor and \code{var}}
\item{penalty}{see \code{\link{lrm}}}
\item{penalty.matrix}{the penalty matrix actually used in the
	estimation}
\item{info.matrix}{a sparse matrix representation of type
	\code{matrix.csr} from the \code{SparseM} package.  This allows the
	full information matrix with all intercepts to be stored efficiently,
	and matrix operations using the Cholesky decomposition to be fast.
	\code{link{vcov.orm}} uses this information to compute the covariance
	matrix for intercepts other than the middle one.}
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Sall J: A monotone regression smoother based on ordinal cumulative
logistic regression, 1991.
	
Le Cessie S, Van Houwelingen JC: Ridge estimators in logistic regression.
Applied Statistics 41:191--201, 1992.

Verweij PJM, Van Houwelingen JC: Penalized likelihood in Cox regression.
Stat in Med 13:2427--2436, 1994.

Gray RJ: Flexible methods for analyzing survival data using splines,
with applications to breast cancer prognosis.  JASA 87:942--951, 1992.

Shao J: Linear model selection by cross-validation.  JASA 88:486--494, 1993.

Verweij PJM, Van Houwelingen JC: Crossvalidation in survival analysis.
Stat in Med 12:2305--2314, 1993.

Harrell FE: Model uncertainty, penalization, and parsimony.  Available
from \url{http://hbiostat.org/talks/iscb98.pdf}.
}
\seealso{
\code{\link{orm.fit}}, \code{\link{predict.orm}}, \code{\link[SparseM:SparseM.solve]{solve}},
\code{\link{rms.trans}}, \code{\link{rms}}, \code{\link[MASS]{polr}},
\code{\link{latex.orm}},  \code{\link{vcov.orm}},
\code{\link[Hmisc]{num.intercepts}},
\code{\link{residuals.orm}}, \code{\link[Hmisc]{na.delete}},
\code{\link[Hmisc]{na.detail.response}},
\code{\link{pentrace}}, \code{\link{rmsMisc}}, \code{\link{vif}},
\code{\link{predab.resample}},
\code{\link{validate.orm}}, \code{\link{calibrate}},
\code{\link{Mean.orm}}, \code{\link{gIndex}}, \code{\link{prModFit}}
}
\examples{
set.seed(1)
n <- 100
y <- round(runif(n), 2)
x1 <- sample(c(-1,0,1), n, TRUE)
x2 <- sample(c(-1,0,1), n, TRUE)
f <- lrm(y ~ x1 + x2, eps=1e-5)
g <- orm(y ~ x1 + x2, eps=1e-5)
max(abs(coef(g) - coef(f)))
w <- vcov(g, intercepts='all') / vcov(f) - 1
max(abs(w))

set.seed(1)
n <- 300
x1 <- c(rep(0,150), rep(1,150))
y <- rnorm(n) + 3*x1
g <- orm(y ~ x1)
g
k <- coef(g)
i <- num.intercepts(g)
h <- orm(y ~ x1, family=probit)
ll <- orm(y ~ x1, family=loglog)
cll <- orm(y ~ x1, family=cloglog)
cau <- orm(y ~ x1, family=cauchit)
x <- 1:i
z <- list(logistic=list(x=x, y=coef(g)[1:i]),
          probit  =list(x=x, y=coef(h)[1:i]),
          loglog  =list(x=x, y=coef(ll)[1:i]),
          cloglog =list(x=x, y=coef(cll)[1:i]))
labcurve(z, pl=TRUE, col=1:4, ylab='Intercept')

tapply(y, x1, mean)
m <- Mean(g)
m(w <- k[1] + k['x1']*c(0,1))
mh <- Mean(h)
wh <- coef(h)[1] + coef(h)['x1']*c(0,1)
mh(wh)

qu <- Quantile(g)
# Compare model estimated and empirical quantiles
cq <- function(y) {
   cat(qu(.1, w), tapply(y, x1, quantile, probs=.1), '\n')
   cat(qu(.5, w), tapply(y, x1, quantile, probs=.5), '\n')
   cat(qu(.9, w), tapply(y, x1, quantile, probs=.9), '\n')
   }
cq(y)

# Try on log-normal model
g <- orm(exp(y) ~ x1)
g
k <- coef(g)
plot(k[1:i])
m <- Mean(g)
m(w <- k[1] + k['x1']*c(0,1))
tapply(exp(y), x1, mean)

qu <- Quantile(g)
cq(exp(y))

# Compare predicted mean with ols for a continuous x
set.seed(3)
n <- 200
x1 <- rnorm(n)
y <- x1 + rnorm(n)
dd <- datadist(x1); options(datadist='dd')
f <- ols(y ~ x1)
g <- orm(y ~ x1, family=probit)
h <- orm(y ~ x1, family=logistic)
w <- orm(y ~ x1, family=cloglog)
mg <- Mean(g); mh <- Mean(h); mw <- Mean(w)
r <- rbind(ols      = Predict(f, conf.int=FALSE),
           probit   = Predict(g, conf.int=FALSE, fun=mg),
           logistic = Predict(h, conf.int=FALSE, fun=mh),
           cloglog  = Predict(w, conf.int=FALSE, fun=mw))
plot(r, groups='.set.')

# Compare predicted 0.8 quantile with quantile regression
qu <- Quantile(g)
qu80 <- function(lp) qu(.8, lp)
f <- Rq(y ~ x1, tau=.8)
r <- rbind(probit   = Predict(g, conf.int=FALSE, fun=qu80),
           quantreg = Predict(f, conf.int=FALSE))
plot(r, groups='.set.')

# Verify transformation invariance of ordinal regression
ga <- orm(exp(y) ~ x1, family=probit)
qua <- Quantile(ga)
qua80 <- function(lp) log(qua(.8, lp))
r <- rbind(logprobit = Predict(ga, conf.int=FALSE, fun=qua80),
           probit    = Predict(g,  conf.int=FALSE, fun=qu80))
plot(r, groups='.set.')

# Try the same with quantile regression.  Need to transform x1
fa <- Rq(exp(y) ~ rcs(x1,5), tau=.8)
r <- rbind(qr    = Predict(f, conf.int=FALSE),
           logqr = Predict(fa, conf.int=FALSE, fun=log))
plot(r, groups='.set.')
options(datadist=NULL)
\dontrun{
## Simulate power and type I error for orm logistic and probit regression
## for likelihood ratio, Wald, and score chi-square tests, and compare
## with t-test
require(rms)
set.seed(5)
nsim <- 2000
r <- NULL
for(beta in c(0, .4)) {
  for(n in c(10, 50, 300)) {
    cat('beta=', beta, '  n=', n, '\n\n')
    plogistic <- pprobit <- plogistics <- pprobits <- plogisticw <-
      pprobitw <- ptt <- numeric(nsim)
    x <- c(rep(0, n/2), rep(1, n/2))
    pb <- setPb(nsim, every=25, label=paste('beta=', beta, '  n=', n))
    for(j in 1:nsim) {
      pb(j)
      y <- beta*x + rnorm(n)
      tt <- t.test(y ~ x)
      ptt[j] <- tt$p.value
      f <- orm(y ~ x)
      plogistic[j]  <- f$stats['P']
      plogistics[j] <- f$stats['Score P']
      plogisticw[j] <- 1 - pchisq(coef(f)['x']^2 / vcov(f)[2,2], 1)
      f <- orm(y ~ x, family=probit)
      pprobit[j]  <- f$stats['P']
      pprobits[j] <- f$stats['Score P']
      pprobitw[j] <- 1 - pchisq(coef(f)['x']^2 / vcov(f)[2,2], 1)
    }
    if(beta == 0) plot(ecdf(plogistic))
    r <- rbind(r, data.frame(beta         = beta, n=n,
                             ttest        = mean(ptt        < 0.05),
                             logisticlr   = mean(plogistic  < 0.05),
                             logisticscore= mean(plogistics < 0.05),
                             logisticwald = mean(plogisticw < 0.05),
                             probit       = mean(pprobit    < 0.05),
                             probitscore  = mean(pprobits   < 0.05),
                             probitwald   = mean(pprobitw   < 0.05)))
  }
}
print(r)
#  beta   n  ttest logisticlr logisticscore logisticwald probit probitscore probitwald
#1  0.0  10 0.0435     0.1060        0.0655        0.043 0.0920      0.0920     0.0820
#2  0.0  50 0.0515     0.0635        0.0615        0.060 0.0620      0.0620     0.0620
#3  0.0 300 0.0595     0.0595        0.0590        0.059 0.0605      0.0605     0.0605
#4  0.4  10 0.0755     0.1595        0.1070        0.074 0.1430      0.1430     0.1285
#5  0.4  50 0.2950     0.2960        0.2935        0.288 0.3120      0.3120     0.3120
#6  0.4 300 0.9240     0.9215        0.9205        0.920 0.9230      0.9230     0.9230
}
}
\keyword{category}
\keyword{models}
\concept{logistic regression model}
\concept{ordinal logistic model}
\concept{proportional odds model}
\concept{ordinal response}
