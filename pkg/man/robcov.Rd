\name{robcov}
\alias{robcov}
\title{Robust Covariance Matrix Estimates}

\description{
Uses the Huber-White method to adjust the variance-covariance matrix of
a fit from maximum likelihood or least squares, to correct for
heteroscedasticity and for correlated responses from cluster samples.
The method uses the ordinary estimates of regression coefficients and
other parameters of the model, but involves correcting the covariance
matrix for model misspecification and sampling design. 
Models currently implemented are models that have a 
\code{residuals(fit,type="score")} function implemented, such as \code{lrm}, 
\code{cph}, \code{coxph}, and ordinary linear models (\code{ols}).
The fit must have specified the \code{x=TRUE} and \code{y=TRUE} options for certain models.
Observations in different clusters are assumed to be independent.
For the special case where every cluster contains one observation, the
corrected covariance matrix returned is the "sandwich" estimator
(see Lin and Wei). This is a consistent estimate of the covariance matrix
even if the model is misspecified (e.g. heteroscedasticity, underdispersion,
wrong covariate form).


For the special case of ols fits, \code{robcov} can compute the improved
(especially for small samples) Efron estimator that adjusts for
natural heterogeneity of residuals (see Long and Ervin (2000)
estimator HC3).
}
\usage{
robcov(fit, cluster, method=c('huber','efron'))
}
\arguments{
\item{fit}{
a fit object from the \code{rms} series
}
\item{cluster}{
a variable indicating groupings. \code{cluster} may be any type of vector
(factor, character, integer).  NAs are not allowed.
Unique values of \code{cluster} indicate
possibly correlated groupings of observations. Note the data used in
the fit and stored in \code{fit$x} and \code{fit$y} may have had observations
containing missing values deleted. It is assumed that if any NAs were
removed during the original model fitting, an \code{naresid} function
exists to restore NAs so that the rows of the score matrix coincide
with \code{cluster}.
If \code{cluster} is omitted,
it defaults to the integers 1,2,\dots,n to obtain the "sandwich" robust
covariance matrix estimate.
}
\item{method}{
can set to \code{"efron"} for ols fits (only).  Default is Huber-White
estimator of the covariance matrix.
}}
\value{
a new fit object with the same class as the original fit,
and with the element \code{orig.var} added. \code{orig.var} is
the covariance matrix of the original fit.  Also, the original \code{var}
component is replaced with the new Huberized estimates.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Huber, PJ. Proc Fifth Berkeley Symposium Math Stat 1:221--33, 1967.


White, H. Econometrica 50:1--25, 1982.


Lin, DY, Wei, LJ. JASA 84:1074--8, 1989.


Rogers, W.  Stata Technical Bulletin STB-8, p. 15--17, 1992.


Rogers, W.  Stata Release 3 Manual, \code{deff}, \code{loneway}, \code{huber}, \code{hreg}, \code{hlogit}
functions.


Long, JS, Ervin, LH.  The American Statistician 54:217--224, 2000.
}
\section{Warnings}{
Adjusted \code{ols} fits do not have the corrected standard errors printed with
\code{print.ols}. Use \code{sqrt(diag(adjfit$var))} to get this, where
\code{adjfit} is the result of \code{robcov}.
}
\seealso{
\code{\link{bootcov}}, \code{\link{naresid}}, \code{\link{residuals.cph}}
}
\examples{
# In OLS test against more manual approach
set.seed(1)
n <- 15
x1 <- 1:n
x2 <- sample(1:n)
y <- round(x1 + x2 + 8*rnorm(n))
f <- ols(y ~ x1 + x2, x=TRUE, y=TRUE)
vcov(f)
vcov(robcov(f))
X <- f$x
G <- diag(resid(f)^2)
solve(t(X) \%*\% X) \%*\% (t(X) \%*\% G \%*\% X) \%*\% solve(t(X) \%*\% X)

# Duplicate data and adjust for intra-cluster correlation to see that
# the cluster sandwich estimator completely ignored the duplicates
x1 <- c(x1,x1)
x2 <- c(x2,x2)
y  <- c(y, y)
g <- ols(y ~ x1 + x2, x=TRUE, y=TRUE)
vcov(robcov(g, c(1:n, 1:n)))

# A dataset contains a variable number of observations per subject,
# and all observations are laid out in separate rows. The responses
# represent whether or not a given segment of the coronary arteries
# is occluded. Segments of arteries may not operate independently
# in the same patient.  We assume a "working independence model" to
# get estimates of the coefficients, i.e., that estimates assuming
# independence are reasonably efficient.  The job is then to get
# unbiased estimates of variances and covariances of these estimates.

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
g <- robcov(f, id)
diag(g$var)/diag(f$var)
# add ,group=w to re-sample from within each level of w
anova(g)            # cluster-adjusted Wald statistics
# fastbw(g)         # cluster-adjusted backward elimination
plot(Predict(g, age=30:70, sex='female'))  # cluster-adjusted confidence bands


# Get design effects based on inflation of the variances when compared
# with bootstrap estimates which ignore clustering
g2 <- robcov(f)
diag(g$var)/diag(g2$var)


# Get design effects based on pooled tests of factors in model
anova(g2)[,1] / anova(g)[,1]




# A dataset contains one observation per subject, but there may be
# heteroscedasticity or other model misspecification. Obtain
# the robust sandwich estimator of the covariance matrix.


# f <- ols(y ~ pol(age,3), x=TRUE, y=TRUE)
# f.adj <- robcov(f)
}
\keyword{models}
\keyword{regression}
\keyword{robust}
\concept{cluster sampling}
\concept{intra-class correlation}
