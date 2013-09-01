\name{Function}
\alias{Function.rms}
\alias{Function.cph}
\alias{sascode}
\alias{perlcode}
\title{Compose an S Function to Compute X beta from a Fit}
\description{
\code{Function} is a class of functions for creating other S functions.
\code{Function.rms} is the method for creating S functions to compute
X beta, based on a model fitted with \code{rms} in effect.  
Like \code{latexrms}, \code{Function.rms} simplifies restricted cubic
spline functions and factors out terms in second-order interactions.
\code{Function.rms} will not work for models that have third-order
interactions involving restricted cubic splines.
\code{Function.cph} is a particular method for handling fits from
\code{cph}, for which an intercept (the negative of the centering
constant) is added to  
the model.  \code{sascode} is a function that takes an S function such
as one created by \code{Function} and does most of the editing
to turn the function definition into
a fragment of SAS code for computing X beta from the fitted model, along
with assignment statements that initialize predictors to reference
values.
\code{perlcode} similarly creates Perl code to evaluate a fitted
regression model.
}
\usage{
\method{Function}{rms}(object, intercept=NULL,
 digits=max(8, .Options$digits), \dots)
\method{Function}{cph}(object, intercept=-object$center, \dots)

# Use result as fun(predictor1=value1, predictor2=value2, \dots)

sascode(object, file='', append=FALSE)

perlcode(object)
}
\arguments{
\item{object}{
a fit created with \code{rms} in effect
}
\item{intercept}{
an intercept value to use (not allowed to be specified to \code{Function.cph}).
The intercept is usually retrieved from the regression coefficients
automatically.
}
\item{digits}{
number of significant digits to use for coefficients and knot locations
}
\item{file}{
name of a file in which to write the SAS code.  Default is to write to
standard output.
}
\item{append}{
set to \code{TRUE} to have \code{sascode} append code to an existing file named
\code{file}.
}
\item{\dots}{arguments to pass to \code{Function.rms} from
  \code{Function.cph}}
}
\value{
\code{Function} returns an S-Plus function that can be invoked in any
usual context.  The function has one argument per predictor variable,
and the default values of the predictors are set to \code{adjust-to} values
(see \code{datadist}).  Multiple predicted X beta values may be calculated
by specifying vectors as arguments to the created function.
All non-scalar argument values must have the same length.
\code{perlcode} returns a character string with embedded newline characters.
}
\author{
Frank Harrell, Jeremy Stephens, and Thomas Dupont\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{latexrms}}, \code{\link[Hmisc]{transcan}},
\code{\link{predict.rms}}, \code{\link{rms}}, \code{\link{rms.trans}} 
}
\examples{
set.seed(1331)
x1 <- exp(rnorm(100))
x2 <- factor(sample(c('a','b'),100,rep=TRUE))
dd <- datadist(x1, x2)
options(datadist='dd')
y  <- log(x1)^2+log(x1)*(x2=='b')+rnorm(100)/4
f  <- ols(y ~ pol(log(x1),2)*x2)
f$coef
g  <- Function(f, digits=5)
g
sascode(g)
cat(perlcode(g), '\n')
g()
g(x1=c(2,3), x2='b')   #could omit x2 since b is default category
predict(f, expand.grid(x1=c(2,3),x2='b'))
g8 <- Function(f)   # default is 8 sig. digits
g8(x1=c(2,3), x2='b')
options(datadist=NULL)


\dontrun{
# Make self-contained functions for computing survival probabilities
# using a log-normal regression
f <- psm(Srv(d.time, death) ~ rcs(age,4)*sex, dist='gaussian')
g <- Function(f)
surv <- Survival(f)
# Compute 2 and 5-year survival estimates for 50 year old male
surv(c(2,5), g(age=50, sex='male'))
}
}
\keyword{regression}
\keyword{methods}
\keyword{interface}
\keyword{models}
\keyword{survival}
\keyword{math}
\concept{logistic regression model}
