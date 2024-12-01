\name{bootcov}
\alias{bootcov}
\alias{bootplot}
\alias{bootplot.bootcov}
\alias{confplot}
\alias{confplot.bootcov}
\alias{histdensity}
\title{Bootstrap Covariance and Distribution for Regression Coefficients}
\description{
\code{bootcov} computes a bootstrap estimate of the covariance matrix for a set
of regression coefficients from \code{ols}, \code{lrm}, \code{cph},
\code{psm}, \code{Rq}, and any
other fit where \code{x=TRUE, y=TRUE} was used to store the data used in making
the original regression fit and where an appropriate \code{fitter} function
is provided here.  The estimates obtained are not conditional on
the design matrix, but are instead unconditional estimates.  For
small sample sizes, this will make a difference as the unconditional
variance estimates are larger.  This function will also obtain
bootstrap estimates corrected for cluster sampling (intra-cluster
correlations) when a "working independence" model was used to fit
data which were correlated within clusters.  This is done by substituting
cluster sampling with replacement for the usual simple sampling with
replacement.  \code{bootcov} has an option (\code{coef.reps}) that causes all
of the regression coefficient estimates from all of the bootstrap
re-samples to be saved, facilitating computation of nonparametric
bootstrap confidence limits and plotting of the distributions of the
coefficient estimates (using histograms and kernel smoothing estimates).


The \code{loglik} option facilitates the calculation of simultaneous
confidence regions from quantities of interest that are functions of
the regression coefficients, using the method of Tibshirani(1996).
With Tibshirani's method, one computes the objective criterion (-2 log
likelihood evaluated at the bootstrap estimate of \eqn{\beta}{beta} but with
respect to the original design matrix and response vector) for the
original fit as well as for all of the bootstrap fits.  The confidence
set of the regression coefficients is the set of all coefficients that
are associated with objective function values that are less than or
equal to say the 0.95 quantile of the vector of \code{B + 1} objective
function values.  For the coefficients satisfying this condition,
predicted values are computed at a user-specified design matrix \code{X},
and minima and maxima of these predicted values (over the qualifying
bootstrap repetitions) are computed to derive the final simultaneous
confidence band.


The \code{bootplot} function takes the output of \code{bootcov} and 
either plots a histogram and kernel density
estimate of specified regression coefficients (or linear combinations
of them through the use of a specified design matrix \code{X}), or a
\code{qqnorm} plot of the quantities of interest to check for normality of
the maximum likelihood estimates.  \code{bootplot} draws vertical lines at
specified quantiles of the bootstrap distribution, and returns these
quantiles for possible printing by the user.  Bootstrap estimates may
optionally be transformed by a user-specified function \code{fun} before
plotting.


The \code{confplot} function also uses the output of \code{bootcov} but to
compute and optionally plot nonparametric bootstrap pointwise confidence
limits or (by default) Tibshirani (1996) simultaneous confidence sets.
A design matrix must be specified to allow \code{confplot} to compute
quantities of interest such as predicted values across a range
of values or differences in predicted values (plots of effects of
changing one or more predictor variable values).


\code{bootplot} and \code{confplot} are actually generic functions, with
the particular functions \code{bootplot.bootcov} and \code{confplot.bootcov}
automatically invoked for \code{bootcov} objects.


A service function called \code{histdensity} is also provided (for use with
\code{bootplot}).  It runs \code{hist} and \code{density} on the same plot, using
twice the number of classes than the default for \code{hist}, and 1.5 times the
\code{width} than the default used by \code{density}.


A comprehensive example demonstrates the use of all of the functions.
}
\usage{
bootcov(fit, cluster, B=200, fitter, 
        coef.reps=TRUE, loglik=FALSE,
        pr=FALSE, maxit=15, eps=0.0001, group=NULL, stat=NULL,
        seed=sample(10000, 1), ...)


bootplot(obj, which=1 : ncol(Coef), X,
         conf.int=c(.9,.95,.99),
         what=c('density', 'qqnorm', 'box'),
         fun=function(x) x, labels., \dots)


confplot(obj, X, against, 
         method=c('simultaneous','pointwise'),
         conf.int=0.95, fun=function(x)x,
         add=FALSE, lty.conf=2, \dots)


histdensity(y, xlab, nclass, width, mult.width=1, \dots)
}
\arguments{
\item{fit}{
a fit object containing components \code{x} and \code{y}.  For fits from
\code{cph}, the \code{"strata"} attribute of the \code{x} component is used to
obtain the vector of stratum codes.
}
\item{obj}{
an object created by \code{bootcov} with \code{coef.reps=TRUE}.
}
\item{X}{
a design matrix specified to \code{confplot}.  See \code{predict.rms} or
\code{contrast.rms}.  For \code{bootplot}, \code{X} is optional.
}
\item{y}{
a vector to pass to \code{histdensity}.  \code{NA}s are ignored.
}
\item{cluster}{
a variable indicating groupings. \code{cluster} may be any type of vector
(factor, character, integer).
Unique values of \code{cluster} indicate
possibly correlated groupings of observations. Note the data used in
the fit and stored in \code{fit$x} and \code{fit$y} may have had observations
containing missing values deleted.  It is assumed that if there were
any NAs, an \code{naresid} function exists for the class of \code{fit}. This
function restores NAs so that the rows of the design matrix
coincide with \code{cluster}.
}
\item{B}{
number of bootstrap repetitions.  Default is 200.
}
\item{fitter}{
the name of a function with arguments \code{(x,y)} that will fit bootstrap
samples.  Default is taken from the class of \code{fit} if it is
\code{ols}, \code{lrm}, \code{cph}, \code{psm}, \code{Rq}.
}
\item{coef.reps}{
set to \code{TRUE} if you want to store a matrix of all bootstrap regression
coefficient estimates in the returned component \code{boot.Coef}.
}
\item{loglik}{
set to \code{TRUE} to store -2 log likelihoods for each bootstrap model,
evaluated against the original \code{x} and \code{y} data.  The default
is to do this when \code{coef.reps} is specified as \code{TRUE}.  The
use of \code{loglik=TRUE} assumes that an \code{oos.loglik} method
exists for the type of model being analyzed, 
to calculate out-of-sample -2 log likelihoods (see \code{rmsMisc}).
After the \code{B} -2 log likelihoods (stored in the element named
\code{boot.loglik} in the returned fit object), the \code{B+1} element is
the -2 log likelihood for the original model fit.
}
\item{pr}{
set to \code{TRUE} to print the current sample number to monitor progress.
}
\item{maxit}{maximum number of iterations, to pass to \code{fitter}}
\item{eps}{argument to pass to various fitters}
\item{group}{
a grouping variable used to stratify the sample upon bootstrapping.
This allows one to handle k-sample problems, i.e., each bootstrap
sample will be forced to select the same number of observations from
each level of group as the number appearing in the original dataset.
You may specify both \code{group} and \code{cluster}.
}
\item{stat}{
  a single character string specifying the name of a \code{stats}
  element produced by the fitting function to save over the bootstrap
  repetitions.  The vector of saved statistics will be in the
  \code{boot.stats} part of the list returned by \code{bootcov}.
}
\item{seed}{random number seed for \code{set.seed}, defaults to a random
	integer between 1 and 10000; user should specify a constant for reproducibility}
\item{which}{
one or more integers specifying which regression coefficients to
plot for \code{bootplot}
}
\item{conf.int}{
a vector (for \code{bootplot}, default is \code{c(.9,.95,.99)}) or scalar 
(for \code{confplot}, default is \code{.95}) confidence level.
}
\item{what}{
for \code{bootplot}, specifies whether a density or a q-q plot is made,
a \code{ggplot2} is used to produce a box plot of all coefficients over
the bootstrap reps
}
\item{fun}{
for \code{bootplot} or \code{confplot} specifies a function used to translate
the quantities of interest before analysis.  A common choice is
\code{fun=exp} to compute anti-logs, e.g., odds ratios.
}
\item{labels.}{
a vector of labels for labeling the axes in plots produced by \code{bootplot}.
Default is row names of \code{X} if there are any, or sequential integers.
}
\item{\dots}{
For \code{bootcov}, extra arguments to pass to \code{lrm.fit} for logistic models.
For \code{bootplot} these are optional arguments passed to
\code{histdensity}.  Also may be optional arguments passed to
\code{plot} by \code{confplot} or optional arguments passed to
\code{hist} from \code{histdensity}, such as \code{xlim} and
\code{breaks}.  The argument \code{probability=TRUE} is always passed to
\code{hist}.
}
\item{against}{
For \code{confplot}, specifying \code{against} causes a plot to be made (or added to).
The \code{against} variable is associated with rows of \code{X} and is used as the
x-coordinates.
}
\item{method}{
specifies whether \code{"pointwise"} or \code{"simultaneous"} confidence regions
are derived by \code{confplot}.  The default is simultaneous.
}
\item{add}{
set to \code{TRUE} to add to an existing plot, for \code{confplot}
}
\item{lty.conf}{
line type for plotting confidence bands in \code{confplot}.  Default is
2 for dotted lines.
}
\item{xlab}{
label for x-axis for \code{histdensity}.  Default is \code{label} attribute or
argument name if there is no \code{label}.
}
\item{nclass}{
passed to \code{hist} if present
}
\item{width}{
passed to \code{density} if present
}
\item{mult.width}{
multiplier by which to adjust the default \code{width} passed to \code{density}.
Default is 1.
}
}
\value{
a new fit object with class of the original object and with the element
\code{orig.var} added. \code{orig.var} is
the covariance matrix of the original fit.  Also, the original \code{var}
component is replaced with the new bootstrap estimates.  The component
\code{boot.coef} is also added.  This contains the mean bootstrap estimates
of regression coefficients (with a log scale element added if
applicable).  \code{boot.Coef} is added if \code{coef.reps=TRUE}.
\code{boot.loglik} is added if \code{loglik=TRUE}.  If \code{stat} is
specified an additional vector \code{boot.stats} will be contained in
the returned object.  \code{B} contains the number of successfully fitted
bootstrap resamples.    A component
\code{clusterInfo} is added to contain elements \code{name} and \code{n}
holding the name of the \code{cluster} variable and the number of clusters.


\code{bootplot} returns a (possible matrix) of quantities of interest and
the requested quantiles of them.  \code{confplot} returns three vectors:
\code{fitted}, \code{lower}, and \code{upper}.
}
\section{Side Effects}{
\code{bootcov} prints if \code{pr=TRUE}
}
\details{
If the fit has a scale parameter (e.g., a fit from \code{psm}), the log
of the individual bootstrap scale estimates are added to the vector
of parameter estimates and and column and row for the log scale are
added to the new covariance matrix (the old covariance matrix also
has this row and column).

For \code{Rq} fits, the \code{tau}, \code{method}, and \code{hs}
arguments are taken from the original fit.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
\email{fh@fharrell.com}\cr

Bill Pikounis\cr
Biometrics Research Department\cr
Merck Research Laboratories\cr
\url{https://billpikounis.com/wpb/}
}
\references{
Feng Z, McLerran D, Grizzle J (1996): A comparison of statistical methods for
clustered data analysis with Gaussian error.  Stat in Med 15:1793--1806.

Tibshirani R, Knight K (1996): Model search and inference by bootstrap 
"bumping". Department of Statistics, University of Toronto.  Technical
report available from
\cr
http://www-stat.stanford.edu/~tibs/.
Presented at the Joint Statistical Meetings,
Chicago, August 1996.
}
\seealso{
\code{\link{robcov}}, \code{\link{sample}}, \code{\link{rms}},
\code{\link{lm.fit}}, \code{\link{lrm.fit}},
\code{\link[survival]{survival-internal}}, 
\code{\link{predab.resample}}, \code{\link{rmsMisc}},
\code{\link{Predict}}, \code{\link{gendata}}, 
\code{\link{contrast.rms}}, \code{\link{Predict}}, \code{\link{setPb}},
\code{multiwayvcov::cluster.boot}
}
\examples{
set.seed(191)
x <- exp(rnorm(200))
logit <- 1 + x/2
y <- ifelse(runif(200) <= plogis(logit), 1, 0)
f <- lrm(y ~ pol(x,2), x=TRUE, y=TRUE)
g <- bootcov(f, B=50, pr=TRUE, seed=3)
anova(g)    # using bootstrap covariance estimates
fastbw(g)   # using bootstrap covariance estimates
beta <- g$boot.Coef[,1]
hist(beta, nclass=15)     #look at normality of parameter estimates
qqnorm(beta)
# bootplot would be better than these last two commands


# A dataset contains a variable number of observations per subject,
# and all observations are laid out in separate rows. The responses
# represent whether or not a given segment of the coronary arteries
# is occluded. Segments of arteries may not operate independently
# in the same patient.  We assume a "working independence model" to
# get estimates of the coefficients, i.e., that estimates assuming
# independence are reasonably efficient.  The job is then to get
# unbiased estimates of variances and covariances of these estimates.


set.seed(2)
n.subjects <- 30
ages <- rnorm(n.subjects, 50, 15)
sexes  <- factor(sample(c('female','male'), n.subjects, TRUE))
logit <- (ages-50)/5
prob <- plogis(logit)  # true prob not related to sex
id <- sample(1:n.subjects, 300, TRUE) # subjects sampled multiple times
table(table(id))  # frequencies of number of obs/subject
age <- ages[id]
sex <- sexes[id]
# In truth, observations within subject are independent:
y   <- ifelse(runif(300) <= prob[id], 1, 0)
f <- lrm(y ~ lsp(age,50)*sex, x=TRUE, y=TRUE)
g <- bootcov(f, id, B=50, seed=3)  # usually do B=200 or more
diag(g$var)/diag(f$var)
# add ,group=w to re-sample from within each level of w
anova(g)            # cluster-adjusted Wald statistics
# fastbw(g)         # cluster-adjusted backward elimination
plot(Predict(g, age=30:70, sex='female'))  # cluster-adjusted confidence bands


# Get design effects based on inflation of the variances when compared
# with bootstrap estimates which ignore clustering
g2 <- bootcov(f, B=50, seed=3)
diag(g$var)/diag(g2$var)


# Get design effects based on pooled tests of factors in model
anova(g2)[,1] / anova(g)[,1]


# Simulate binary data where there is a strong 
# age x sex interaction with linear age effects 
# for both sexes, but where not knowing that
# we fit a quadratic model.  Use the bootstrap
# to get bootstrap distributions of various
# effects, and to get pointwise and simultaneous
# confidence limits


set.seed(71)
n   <- 500
age <- rnorm(n, 50, 10)
sex <- factor(sample(c('female','male'), n, rep=TRUE))
L   <- ifelse(sex=='male', 0, .1*(age-50))
y   <- ifelse(runif(n)<=plogis(L), 1, 0)


f <- lrm(y ~ sex*pol(age,2), x=TRUE, y=TRUE)
b <- bootcov(f, B=50, loglik=TRUE, pr=TRUE, seed=3)   # better: B=500


par(mfrow=c(2,3))
# Assess normality of regression estimates
bootplot(b, which=1:6, what='qq')
# They appear somewhat non-normal


# Plot histograms and estimated densities 
# for 6 coefficients
w <- bootplot(b, which=1:6)
# Print bootstrap quantiles
w$quantiles

# Show box plots for bootstrap reps for all coefficients
bootplot(b, what='box')


# Estimate regression function for females
# for a sequence of ages
ages <- seq(25, 75, length=100)
label(ages) <- 'Age'


# Plot fitted function and pointwise normal-
# theory confidence bands
par(mfrow=c(1,1))
p <- Predict(f, age=ages, sex='female')
plot(p)
# Save curve coordinates for later automatic
# labeling using labcurve in the Hmisc library
curves <- vector('list',8)
curves[[1]] <- with(p, list(x=age, y=lower))
curves[[2]] <- with(p, list(x=age, y=upper))


# Add pointwise normal-distribution confidence 
# bands using unconditional variance-covariance
# matrix from the 500 bootstrap reps
p <- Predict(b, age=ages, sex='female')
curves[[3]] <- with(p, list(x=age, y=lower))
curves[[4]] <- with(p, list(x=age, y=upper))


dframe <- expand.grid(sex='female', age=ages)
X <- predict(f, dframe, type='x')  # Full design matrix


# Add pointwise bootstrap nonparametric 
# confidence limits
p <- confplot(b, X=X, against=ages, method='pointwise',
              add=TRUE, lty.conf=4)
curves[[5]] <- list(x=ages, y=p$lower)
curves[[6]] <- list(x=ages, y=p$upper)


# Add simultaneous bootstrap confidence band
p <- confplot(b, X=X, against=ages, add=TRUE, lty.conf=5)
curves[[7]] <- list(x=ages, y=p$lower)
curves[[8]] <- list(x=ages, y=p$upper)
lab <- c('a','a','b','b','c','c','d','d')
labcurve(curves, lab, pl=TRUE)


# Now get bootstrap simultaneous confidence set for
# female:male odds ratios for a variety of ages


dframe <- expand.grid(age=ages, sex=c('female','male'))
X <- predict(f, dframe, type='x')  # design matrix
f.minus.m <- X[1:100,] - X[101:200,]
# First 100 rows are for females.  By subtracting
# design matrices are able to get Xf*Beta - Xm*Beta
# = (Xf - Xm)*Beta


confplot(b, X=f.minus.m, against=ages,
         method='pointwise', ylab='F:M Log Odds Ratio')
confplot(b, X=f.minus.m, against=ages,
         lty.conf=3, add=TRUE)


# contrast.rms makes it easier to compute the design matrix for use
# in bootstrapping contrasts:


f.minus.m <- contrast(f, list(sex='female',age=ages),
                         list(sex='male',  age=ages))$X
confplot(b, X=f.minus.m)


# For a quadratic binary logistic regression model use bootstrap
# bumping to estimate coefficients under a monotonicity constraint
set.seed(177)
n <- 400
x <- runif(n)
logit <- 3*(x^2-1)
y <- rbinom(n, size=1, prob=plogis(logit))
f <- lrm(y ~ pol(x,2), x=TRUE, y=TRUE)
k <- coef(f)
k
vertex <- -k[2]/(2*k[3])
vertex


# Outside [0,1] so fit satisfies monotonicity constraint within
# x in [0,1], i.e., original fit is the constrained MLE


g <- bootcov(f, B=50, coef.reps=TRUE, loglik=TRUE, seed=3)
bootcoef <- g$boot.Coef    # 100x3 matrix
vertex <- -bootcoef[,2]/(2*bootcoef[,3])
table(cut2(vertex, c(0,1)))
mono <- !(vertex >= 0 & vertex <= 1)
mean(mono)    # estimate of Prob{monotonicity in [0,1]}


var(bootcoef)   # var-cov matrix for unconstrained estimates
var(bootcoef[mono,])   # for constrained estimates


# Find second-best vector of coefficient estimates, i.e., best
# from among bootstrap estimates
g$boot.Coef[order(g$boot.loglik[-length(g$boot.loglik)])[1],]
# Note closeness to MLE

\dontrun{
# Get the bootstrap distribution of the difference in two ROC areas for
# two binary logistic models fitted on the same dataset.  This analysis
# does not adjust for the bias ROC area (C-index) due to overfitting.
# The same random number seed is used in two runs to enforce pairing.

set.seed(17)
x1 <- rnorm(100)
x2 <- rnorm(100)
y <- sample(0:1, 100, TRUE)
f <- lrm(y ~ x1, x=TRUE, y=TRUE)
g <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
f <- bootcov(f, stat='C', seed=4)
g <- bootcov(g, stat='C', seed=4)
dif <- g$boot.stats - f$boot.stats
hist(dif)
quantile(dif, c(.025,.25,.5,.75,.975))
# Compute a z-test statistic.  Note that comparing ROC areas is far less
# powerful than likelihood or Brier score-based methods
z <- (g$stats['C'] - f$stats['C'])/sd(dif)
names(z) <- NULL
c(z=z, P=2*pnorm(-abs(z)))
}
}
\keyword{models}
\keyword{regression}
\keyword{htest}
\keyword{methods}
\keyword{hplot}
\concept{bootstrap}
\concept{sampling}
