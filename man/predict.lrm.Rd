\name{predict.lrm}
\alias{predict.lrm}
\alias{predict.orm}
\alias{Mean.lrm}
\alias{Mean.orm}
\title{
Predicted Values for Binary and Ordinal Logistic Models
}
\description{
Computes a variety of types of predicted values for fits from
\code{lrm} and \code{orm}, either from the original dataset or for new
observations.  The \code{Mean.lrm} and \code{Mean.orm} functions produce
an R function to compute the predicted mean of a numeric ordered
response variable given the linear predictor, which is assumed to use
the first intercept when it was computed.  The returned function has two
optional arguments if confidence intervals are desired: \code{conf.int}
and the design matrix \code{X}.  When this derived function is called
with nonzero \code{conf.int}, an attribute named \code{limits} is attached
	to the estimated mean.  This is a list with elements \code{lower} and
	\code{upper} containing normal approximations for confidence limits
	using the delta method.
}
\usage{
\method{predict}{lrm}(object, \dots, type=c("lp", "fitted",
            "fitted.ind", "mean", "x", "data.frame",
            "terms", "cterms", "ccterms", "adjto","adjto.data.frame", 
            "model.frame"), se.fit=FALSE, codes=FALSE)
\method{predict}{orm}(object, \dots, type=c("lp", "fitted",
            "fitted.ind", "mean", "x", "data.frame",
            "terms", "cterms", "ccterms", "adjto","adjto.data.frame", 
            "model.frame"), se.fit=FALSE, codes=FALSE)

\method{Mean}{lrm}(object, codes=FALSE, \dots)
\method{Mean}{orm}(object, codes=FALSE, \dots)
}
\arguments{
\item{object}{a object created by \code{lrm} or \code{orm}}
\item{\dots}{
arguments passed to \code{predictrms}, such as \code{kint} and \code{newdata}
(which is used if you are predicting \code{out of data}).  See
\code{predictrms} to see how NAs are handled.  Ignored for other functions.
}
\item{type}{
See \code{predict.rms} for \code{"x", "data.frame", "terms", "cterms",
  "ccterms", "adjto", "adjto.data.frame"} and
\code{"model.frame"}. \code{type="lp"} is used to get 
linear predictors (using the first intercept by default; specify
\code{kint} to use others). \code{type="fitted"}
is used to get all the probabilities \eqn{Y\geq	j}.
\code{type="fitted.ind"} gets all the individual probabilities 
\eqn{Y=j} (not recommended for \code{orm} fits). For an ordinal response
variable, \code{type="mean"} computes 
the estimated mean \eqn{Y} by summing values of \eqn{Y} 
multiplied by the estimated \eqn{Prob(Y=j)}. If \eqn{Y} was a character or
\code{factor} object, the levels are the character values or factor levels,
so these must be translatable to numeric, unless \code{codes=TRUE}.
See the Hannah and Quigley reference below for the method of estimating
(and presenting) the mean score.  If you specify
\code{type="fitted","fitted.ind","mean"} you may not specify \code{kint}.
}
\item{se.fit}{
applies only to \code{type="lp"}, to get standard errors.
}
\item{codes}{
if \code{TRUE}, \code{type="mean"}, \code{Mean.lrm}, and \code{Mean.orm}
use the integer codes \eqn{1,2,\ldots,k} for the \eqn{k}-level response
in computing the predicted mean response.
}
}
\value{
a vector (\code{type="lp"} with \code{se.fit=FALSE}, or
\code{type="mean"} or only one 
observation being predicted), a list (with elements \code{linear.predictors}
and \code{se.fit} if \code{se.fit=TRUE}), a matrix (\code{type="fitted"}
or \code{type="fitted.ind"}), a data frame, or a design matrix.  For
\code{Mean.lrm} and \code{Mean.orm}, the result is an R function.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com\cr
For the \code{Quantile} function:\cr
Qi Liu and Shengxin Tu\cr
Department of Biostatistics, Vanderbilt University
}
\references{
Hannah M, Quigley P: Presentation of ordinal regression analysis on the
original scale.  Biometrics 52:771--5; 1996.
}
\seealso{
\code{\link{lrm}}, \code{\link{orm}}, \code{\link{predict.rms}},
\code{\link{naresid}}, \code{\link{contrast.rms}} 
}
\examples{
# See help for predict.rms for several binary logistic
# regression examples


# Examples of predictions from ordinal models
set.seed(1)
y <- factor(sample(1:3, 400, TRUE), 1:3, c('good','better','best'))
x1 <- runif(400)
x2 <- runif(400)
f <- lrm(y ~ rcs(x1,4)*x2, x=TRUE)     #x=TRUE needed for se.fit
# Get 0.95 confidence limits for Prob[better or best]
L <- predict(f, se.fit=TRUE)           #omitted kint= so use 1st intercept
plogis(with(L, linear.predictors + 1.96*cbind(-se.fit,se.fit)))
predict(f, type="fitted.ind")[1:10,]   #gets Prob(better) and all others
d <- data.frame(x1=c(.1,.5),x2=c(.5,.15))
predict(f, d, type="fitted")        # Prob(Y>=j) for new observation
predict(f, d, type="fitted.ind")    # Prob(Y=j)
predict(f, d, type='mean', codes=TRUE) # predicts mean(y) using codes 1,2,3
m <- Mean(f, codes=TRUE)
lp <- predict(f, d)
m(lp)
# Can use function m as an argument to Predict or nomogram to
# get predicted means instead of log odds or probabilities
dd <- datadist(x1,x2); options(datadist='dd')
m
plot(Predict(f, x1, fun=m), ylab='Predicted Mean')
# Note: Run f through bootcov with coef.reps=TRUE to get proper confidence
# limits for predicted means from the prop. odds model
options(datadist=NULL)
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
