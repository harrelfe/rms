\name{ols}
\alias{ols}
\title{Linear Model Estimation Using Ordinary Least Squares}

\description{
Fits the usual weighted or unweighted linear regression model using the
same fitting routines used by \code{lm}, but also storing the variance-covariance
matrix \code{var} and using traditional dummy-variable coding for categorical
factors.  
Also fits unweighted models using penalized least squares, with the same
penalization options as in the \code{lrm} function.  For penalized estimation,
there is a fitter function call \code{lm.pfit}.
}
\usage{
ols(formula, data=environment(formula), weights, subset, na.action=na.delete, 
    method="qr", model=FALSE,
    x=FALSE, y=FALSE, se.fit=FALSE, linear.predictors=TRUE,
    penalty=0, penalty.matrix, tol=.Machine$double.eps, sigma,
    var.penalty=c('simple','sandwich'), \dots)
}
\arguments{
\item{formula}{
an S formula object, e.g. 
\cr
    Y ~ rcs(x1,5)*lsp(x2,c(10,20))
}
\item{data}{
name of an S data frame containing all needed variables.  Omit this to use a
data frame already in the S ``search list''.
}
\item{weights}{an optional vector of weights to be used in the fitting
          process. If specified, weighted least squares is used with
          weights \code{weights} (that is, minimizing \eqn{sum(w*e^2)});
          otherwise ordinary least squares is used.}
\item{subset}{
an expression defining a subset of the observations to use in the fit.  The default
is to use all observations.  Specify for example \code{age>50 & sex="male"} or
\code{c(1:100,200:300)}
respectively to use the observations satisfying a logical expression or those having
row numbers in the given vector.
}
\item{na.action}{
specifies an S function to handle missing data.  The default is the function \code{na.delete},
which causes observations with any variable missing to be deleted.  The main difference
between \code{na.delete} and the S-supplied function \code{na.omit} is that 
\code{na.delete} makes a list
of the number of observations that are missing on each variable in the model.
The \code{na.action} is usally specified by e.g. \code{options(na.action="na.delete")}.
}
\item{method}{
specifies a particular fitting method, or \code{"model.frame"} instead to return the model frame
of the predictor and response variables satisfying any subset or missing value
checks.
}
\item{model}{
default is \code{FALSE}.  Set to \code{TRUE} to return the model frame
as element \code{model} of the fit object.
}
\item{x}{
default is \code{FALSE}.  Set to \code{TRUE} to return the expanded design matrix as element \code{x}
(without intercept indicators) of the
returned fit object.  Set both \code{x=TRUE} if you are going to use
the \code{residuals} function later to return anything other than ordinary residuals.
}
\item{y}{
default is \code{FALSE}.  Set to \code{TRUE} to return the vector of response values 
as element \code{y} of the fit.
}
\item{se.fit}{
default is \code{FALSE}.  Set to \code{TRUE} to compute the estimated standard errors of
the estimate of \eqn{X\beta}{X beta} and store them in element \code{se.fit}
of the fit. 
}
\item{linear.predictors}{
set to \code{FALSE} to cause predicted values not to be stored
}
\item{penalty}{see \code{lrm}}
\item{penalty.matrix}{see \code{lrm}}
\item{tol}{tolerance for information matrix singularity}
\item{sigma}{
If \code{sigma} is given, it is taken as the actual root mean squared error parameter for the model.  Otherwise \code{sigma} is estimated from the data using the usual formulas (except for penalized models).  It is often convenient to specify
\code{sigma=1} for models with no error, when using \code{fastbw} to find an
approximate model that predicts predicted values from the full model with
a given accuracy.
}
\item{var.penalty}{
the type of variance-covariance matrix to be stored in the \code{var}
component of the fit when penalization is used.  The default is the
inverse of the penalized information matrix.  Specify
\code{var.penalty="sandwich"} to use the sandwich estimator (see below
under \code{var}), which limited simulation studies have shown yields
variances estimates that are too low.
}
\item{\dots}{arguments to pass to \code{\link{lm.wfit}} or
  \code{\link{lm.fit}}}
}
\value{
the same objects returned from \code{lm} (unless \code{penalty} or \code{penalty.matrix}
are given - then an
abbreviated list is returned since \code{lm.pfit} is used as a fitter)
plus the design attributes
(see \code{rms}).
Predicted values are always returned, in the element \code{linear.predictors}.
The vectors or matrix stored if \code{y=TRUE} or \code{x=TRUE} have rows deleted according to \code{subset} and
to missing data, and have names or row names that come from the
data frame used as input data.  If \code{penalty} or \code{penalty.matrix} is given, 
the \code{var} matrix
returned is an improved variance-covariance matrix
for the penalized regression coefficient estimates.  If
\code{var.penalty="sandwich"} (not the default, as limited simulation
studies have found it provides variance estimates that are too low) it
is defined as 
\eqn{\sigma^{2} (X'X + P)^{-1} X'X (X'X + P)^{-1}}, where \eqn{P} is 
\code{penalty factors * penalty.matrix}, with a column and row of zeros
added for the
intercept.  When \code{var.penalty="simple"} (the default), \code{var} is
\eqn{\sigma^{2} (X'X + P)^{-1}}.
The returned list has a vector \code{stats} with named elements
\code{n, Model L.R., d.f., R2, g, Sigma}.  \code{Model L.R.} is the model
likelihood ratio \eqn{\chi^2}{chi-square} statistic, and \code{R2} is
\eqn{R^2}.  For penalized estimation, \code{d.f.} is the 
effective degrees of freedom, which is the sum of the elements of another
vector returned, \code{effective.df.diagonal}, minus one for the
intercept.
\code{g} is the \eqn{g}-index.
\code{Sigma} is the penalized maximum likelihood estimate (see below).
}
\details{
For penalized estimation, the penalty factor on the log likelihood is
\eqn{-0.5 \beta' P \beta / \sigma^2}, where \eqn{P} is defined above.
The penalized maximum likelihood estimate (penalized least squares
or ridge estimate) of \eqn{\beta}{beta} is \eqn{(X'X + P)^{-1} X'Y}.
The maximum likelihood estimate of \eqn{\sigma^2} is \eqn{(sse + \beta'
  P \beta) / n}, where
\code{sse} is the sum of squared errors (residuals).
The \code{effective.df.diagonal} vector is the
diagonal of the matrix \eqn{X'X/(sse/n) \sigma^{2} (X'X + P)^{-1}}.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{anova.rms}},
\code{\link{summary.rms}}, \code{\link{predict.rms}},
\code{\link{fastbw}}, \code{\link{validate}}, \code{\link{calibrate}},
\code{\link{Predict}}, \code{\link{specs.rms}}, \code{\link{cph}},
\code{\link{lrm}}, \code{\link{which.influence}}, \code{\link{lm}},
\code{\link{summary.lm}}, \code{\link{print.ols}},
\code{\link{residuals.ols}}, \code{\link{latex.ols}},
\code{\link[Hmisc]{na.delete}}, \code{\link[Hmisc]{na.detail.response}},
\code{\link{datadist}}, \code{\link{pentrace}}, \code{\link{vif}},
\code{\link[Hmisc]{abs.error.pred}}
}
\examples{
set.seed(1)
x1 <- runif(200)
x2 <- sample(0:3, 200, TRUE)
distance <- (x1 + x2/3 + rnorm(200))^2
d <- datadist(x1,x2)
options(datadist="d")   # No d -> no summary, plot without giving all details


f <- ols(sqrt(distance) ~ rcs(x1,4) + scored(x2), x=TRUE)
# could use d <- datadist(f); options(datadist="d") at this point,
# but predictor summaries would not be stored in the fit object for
# use with Predict, summary.rms.  In that case, the original
# dataset or d would need to be accessed later, or all variable values
# would have to be specified to summary, plot
anova(f)
which.influence(f)
summary(f)
summary.lm(f)    # will only work if penalty and penalty.matrix not used


# Fit a complex model and approximate it with a simple one
x1 <- runif(200)
x2 <- runif(200)
x3 <- runif(200)
x4 <- runif(200)
y <- x1 + x2 + rnorm(200)
f    <- ols(y ~ rcs(x1,4) + x2 + x3 + x4)
pred <- fitted(f)   # or predict(f) or f$linear.predictors
f2   <- ols(pred ~ rcs(x1,4) + x2 + x3 + x4, sigma=1)
# sigma=1 prevents numerical problems resulting from R2=1
fastbw(f2, aics=100000)
# This will find the best 1-variable model, best 2-variable model, etc.
# in predicting the predicted values from the original model
options(datadist=NULL)
}
\keyword{models}
\keyword{regression}
