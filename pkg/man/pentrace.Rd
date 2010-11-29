\name{pentrace}
\alias{pentrace}
\alias{plot.pentrace}
\alias{print.pentrace}
\alias{effective.df}
\title{
Trace AIC and BIC vs. Penalty
}
\description{
For an ordinary unpenalized fit from \code{lrm} or \code{ols} and for a vector or list of penalties, 
fits a series of logistic or linear models using penalized maximum likelihood
estimation, and saves the effective degrees of freedom, Akaike Information
Criterion (\eqn{AIC}), Schwarz Bayesian Information Criterion (\eqn{BIC}), and
Hurvich and Tsai's corrected \eqn{AIC} (\eqn{AIC_c}).  Optionally
\code{pentrace} can 
use the \code{nlminb} function to solve for the optimum penalty factor or
combination of factors penalizing different kinds of terms in the model.
The \code{effective.df} function prints the original and effective
degrees of freedom for a penalized fit or for an unpenalized fit and
the best penalization determined from a previous invocation of
\code{pentrace} if \code{method="grid"} (the default).
The effective d.f. is computed separately for each class of terms in
the model (e.g., interaction, nonlinear).
A \code{plot} method exists to plot the results, and a \code{print} method exists
to print the most pertinent components.  Both \eqn{AIC} and \eqn{BIC}
may be plotted if 
there is only one penalty factor type specified in \code{penalty}.  Otherwise,
the first two types of penalty factors are plotted, showing only the \eqn{AIC}.
}
\usage{
pentrace(fit, penalty, penalty.matrix, 
         method=c('grid','optimize'),
         which=c('aic.c','aic','bic'), target.df,
         fitter, pr=FALSE, tol=1e-7,
         keep.coef=FALSE, complex.more=TRUE, verbose=FALSE, maxit=12, subset)

effective.df(fit, object)

\method{print}{pentrace}(x, \dots)

\method{plot}{pentrace}(x, method=c('points','image'), 
     which=c('effective.df','aic','aic.c','bic'), pch=2, add=FALSE, 
     ylim, \dots)
}
\arguments{
\item{fit}{
a result from \code{lrm} or \code{ols} with \code{x=TRUE, y=TRUE} and without using \code{penalty} or
\code{penalty.matrix}
(or optionally using penalization in the case of \code{effective.df})
}
\item{penalty}{
can be a vector or a list.  If it is a vector, all types of terms in
the model will be penalized by the same amount, specified by elements in
\code{penalty}, with a penalty of zero automatically added.  \code{penalty} can
also be a list in the format documented in the \code{lrm} function, except that
elements of the list can be vectors.  The \code{expand.grid} function is
invoked by \code{pentrace} to generate all possible combinations of
penalties.  For example, specifying 
\code{penalty=list(simple=1:2, nonlinear=1:3)} will generate 6 combinations
to try, so that the analyst can attempt to determine whether penalizing
more complex terms in the model more than the linear or categorical
variable terms will be beneficial.  If \code{complex.more=TRUE}, it is assumed
that the variables given in \code{penalty} are listed in order from less
complex to more complex.  With \code{method="optimize"} \code{penalty} specifies
an initial guess for the penalty or penalties.  If all term types are
to be equally penalized, \code{penalty} should be a single number,
otherwise it should be a list containing single numbers as elements,
e.g., \code{penalty=list(simple=1, nonlinear=2)}.  Experience has shown that the optimization algorithm is more likely to find a reasonable solution when the starting value specified in \code{penalty} is too large rather than too small.
}
\item{object}{
an object returned by \code{pentrace}.  For \code{effective.df}, \code{object} can be
omitted if the \code{fit} was penalized.
}
\item{penalty.matrix}{
see \code{lrm}
}
\item{method}{
The default is \code{method="grid"} to print various indexes for all
combinations of penalty parameters given by the user.  Specify
\code{method="optimize"} to have \code{pentrace} use \code{nlminb} to solve for the
combination of penalty parameters that gives the maximum value of the
objective named in \code{which}, or, if \code{target.df} is given, to find the
combination that yields \code{target.df} effective total degrees of freedom
for the model.  When \code{target.df} is specified, \code{method} is set to
\code{"optimize"} automatically.
For \code{plot.pentrace} this parameter applies only if more than one
penalty term-type was used.  The default is to use open triangles
whose sizes are proportional to the ranks of the AICs, plotting the
first two penalty factors respectively on the x and y  axes.  Use
\code{method="image"} to plot an image plot. 
}
\item{which}{
the objective to maximize for either \code{method}.  Default is \code{"aic.c"} (corrected
AIC).
For \code{plot.pentrace}, \code{which} is a vector of names of criteria to show;
default is to plot all 4 types, with effective d.f. in its own separate plot
}
\item{target.df}{
applies only to \code{method="optimize"}.  See \code{method}.  \code{target.df} makes
sense mainly when a single type of penalty factor is specified.
}
\item{fitter}{
a fitting function.  Default is \code{lrm.fit} (\code{lm.pfit} is always used for \code{ols}).
}
\item{pr}{
set to \code{TRUE} to print intermediate results
}
\item{tol}{
tolerance for declaring a matrix singular (see \code{lrm.fit, solvet})
}
\item{keep.coef}{
set to \code{TRUE} to store matrix of regression  coefficients for all the fits (corresponding
to increasing values of \code{penalty}) in object \code{Coefficients} in the
returned list.  Rows correspond to penalties, columns to regression
parameters.
}
\item{complex.more}{
By default if \code{penalty} is a list, combinations of penalties for which
complex terms are penalized less than less complex terms will be
dropped after \code{expand.grid} is invoked.  Set \code{complex.more=FALSE} to
allow more complex terms to be penalized less.  Currently this option
is ignored for \code{method="optimize"}.
}
\item{verbose}{set to \code{TRUE} to print number of intercepts and sum
  of effective degrees of freedom}
\item{maxit}{
maximum number of iterations to allow in a model fit (default=12).
This is passed to the appropriate fitter function with the correct
argument name.  Increase \code{maxit} if you had to when fitting the
original unpenalized model.
}
\item{subset}{
a logical or integer vector specifying rows of the design and response
matrices to subset in fitting models.  This is most useful for
bootstrapping \code{pentrace} to see if the best penalty can be estimated
with little error so that variation due to selecting the optimal
penalty can
be safely ignored when bootstrapping standard errors of regression
coefficients and measures of predictive accuracy.  See an example below.
}
\item{x}{a result from \code{pentrace}}
\item{pch}{used for \code{method="points"}}
\item{add}{
set to \code{TRUE} to add to an existing plot.  In that case, the effective
d.f. plot is not re-drawn, but the AIC/BIC plot is added to.
}
\item{ylim}{
2-vector of y-axis limits for plots other than effective d.f.
}
\item{...}{
other arguments passed to \code{plot}, \code{lines}, or \code{image}
}}
\value{
a list of class \code{"pentrace"}
with elements \code{penalty, df, objective, fit, var.adj, diag, results.all}, and
optionally \code{Coefficients}.
The first 6 elements correspond to the fit that had the best objective
as named in the \code{which} argument, from the sequence of fits tried.
Here \code{fit} is the fit object from \code{fitter} which was a penalized fit,
\code{diag} is the diagonal of the matrix used to compute the effective
d.f., and \code{var.adj} is Gray (1992) Equation 2.9, which is an improved
covariance matrix for the penalized beta. \code{results.all} is a data
frame whose first few variables are the components of \code{penalty} and
whose other columns are \code{df, aic, bic, aic.c}.  \code{results.all} thus
contains a summary of results for all fits attempted.  When
\code{method="optimize"}, only two components are returned: \code{penalty} and
\code{objective}, and the object does not have a class.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Gray RJ: Flexible methods for analyzing survival data using splines,
with applications to breast cancer prognosis.  JASA 87:942--951, 1992.


Hurvich CM, Tsai, CL: Regression and time series model selection in small
samples.  Biometrika 76:297--307, 1989.
}
\seealso{
\code{\link{lrm}}, \code{\link{ols}}, \code{\link[Hmisc]{solvet}}, \code{\link{rmsMisc}}, \code{\link{image}}
}
\examples{
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)


f <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
         x=TRUE, y=TRUE)
p <- pentrace(f, seq(.2,1,by=.05))
plot(p)
p$diag      # may learn something about fractional effective d.f. 
            # for each original parameter
pentrace(f, list(simple=c(0,.2,.4), nonlinear=c(0,.2,.4,.8,1)))


# Bootstrap pentrace 5 times, making a plot of corrected AIC plot with 5 reps
n <- nrow(f$x)
plot(pentrace(f, seq(.2,1,by=.05)), which='aic.c', 
     col=1, ylim=c(30,120)) #original in black
for(j in 1:5)
  plot(pentrace(f, seq(.2,1,by=.05), subset=sample(n,n,TRUE)), 
       which='aic.c', col=j+1, add=TRUE)


# Find penalty giving optimum corrected AIC.  Initial guess is 1.0
# Not implemented yet
# pentrace(f, 1, method='optimize')


# Find penalty reducing total regression d.f. effectively to 5
# pentrace(f, 1, target.df=5)


# Re-fit with penalty giving best aic.c without differential penalization
f <- update(f, penalty=p$penalty)
effective.df(f)
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
\concept{penalized MLE}
\concept{ridge regression}
\concept{shrinkage}
