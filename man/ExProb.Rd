\name{ExProb}
\alias{ExProb}
\alias{ExProb.orm}
\alias{Survival.orm}
\alias{plot.ExProb}
\title{Function Generators For Exceedance and Survival Probabilities}
\description{
	For an \code{orm} object \code{ExProb} generates a function for computing the
	estimates of the function Prob(Y>=y) given one or more values of the
	linear predictor using the reference (median) intercept.  This
	function can optionally be evaluated at only a set of user-specified
	\code{y} values, otherwise a right-step function is returned.  There
	is a plot method for plotting the step functions, and if more than one
	linear predictor was evaluated multiple step functions are drawn.
	\code{ExProb} is especially useful for \code{\link{nomogram}}.  \code{Survival}
  generates a similar function but for computing survival probabilities Prob(Y>y)
  and adding an origin of zero.  Plotting of survival curves is done with a 
  \code{survplot} method.  \code{Survival.orm} merely calls \code{ExProb.orm}
  with argument \code{survival=TRUE}.  For survival estimation when interval
  censoring is present, times are taken as interval midpoints with intervals
  corresponding to intercepts in the model.
	
	Optionally a normal approximation (normality for the linear predictor) for a confidence
	interval for exceedance probabilities will be computed, if
	\code{conf.int > 0} is specified to the function generated from calling
	\code{ExProb} or \code{Survival}.  For \code{ExProb}, a \code{"lims"} attribute is included
	in the result computed by the derived cumulative probability function.  For \code{Survival},
  the result is a data frame if \code{conf.int} is specified or both time and the requested
  linear predictor are varying.  In the data frame the limits are variables \code{lower} and \code{upper}.
}
\usage{
ExProb(object, \dots)

\method{ExProb}{orm}(object, codes = FALSE, ...)

\method{plot}{ExProb}(x, \dots, data=NULL,
                      xlim=NULL, xlab=x$yname, ylab=expression(Prob(Y>=y)),
                      col=par('col'), col.vert='gray85', pch=20,
                      pch.data=21, lwd=par('lwd'), lwd.data=lwd,
                      lty.data=2, key=TRUE)

\method{Survival}{orm}(object, \dots)
}
\arguments{
  \item{object}{a fit object from \code{orm}.  For \code{Survival} the fit may be from \code{orm.fit}.
  This is used to estimate survival curves when there are no predictors in the model.  In the case
  the link function (\code{family} argument to \code{orm.fit}) does not affect survival probabilities
  but does affect confidence limits.  To get the same confidence intervals as \code{survival:survfit.formula} 
  use \code{ormfit(y=, family='loglog')} to correspond to \code{survfit(..., conf.type='log-log')}.}
  \item{codes}{if \code{TRUE}, \code{ExProb} use the integer codes
     \eqn{1,2,\ldots,k} for the \eqn{k}-level response instead of its
     original unique values} 
  \item{\dots}{ignored for \code{ExProb}.  Passed to \code{plot} for
    \code{plot.ExProb}}
	\item{data}{Specify \code{data} if you want to add stratified empirical
    probabilities to the graph.  If \code{data} is a numeric vector, it
    is assumed that no groups are present.  Otherwise \code{data} must
    be a list or data frame where the first variable is the grouping
    variable (corresponding to what made the linear predictor vary) and
		the second variable is the data vector for the \code{y} variable.
		The rows of data should be sorted to be in order of the linear
    predictor argument.
	}
	\item{x}{an object created by running the function created by \code{ExProb}}
  \item{xlim}{limits for x-axis; default is range of observed \code{y}}
  \item{xlab}{x-axis label}
	\item{ylab}{y-axis label}
  \item{col}{color for horizontal lines and points}
	\item{col.vert}{color for vertical discontinuities}
	\item{pch}{plotting symbol for predicted curves}
	\item{lwd}{line width for predicted curves}
	\item{pch.data,lwd.data,lty.data}{plotting parameters for data}
	\item{key}{set to \code{FALSE} to suppress key in plot if \code{data}
                      is given}
}
\value{
	\code{ExProb} and \code{Survival} return an R function.  Running the function returns an
  object of class \code{"ExProb"} for \code{ExProb}, or a data frame or vector for \code{Survival}.
}
\author{Frank Harrell and Shengxin Tu}
\seealso{\code{\link{orm}}, \code{\link{Quantile.orm}}}
\examples{
set.seed(1)
x1 <- runif(200)
yvar <- x1 + runif(200)
f <- orm(yvar ~ x1)
d <- ExProb(f)
lp <- predict(f, newdata=data.frame(x1=c(.2,.8)))
w <- d(lp)
s1 <- abs(x1 - .2) < .1
s2 <- abs(x1 - .8) < .1
plot(w, data=data.frame(x1=c(rep(.2, sum(s1)), rep(.8, sum(s2))),
                        yvar=c(yvar[s1], yvar[s2])))

qu <- Quantile(f)
abline(h=c(.1,.5), col='gray80')
abline(v=qu(.5, lp), col='gray80')
abline(v=qu(.9, lp), col='green')
\dontrun{
  Y <- Ocens(dtime, ifelse(censored, Inf, dtime))
  f <- orm(Y ~ x, family='loglog')
  s <- Survival(f)
  s()  # all times
  s(times=c(1, 3))
  d <- data.frame(x=2:4)
  s(X=predict(f, d, conf.int=0.95)  # all times
  s(lp=predict(f, d))  # same surv estimates, no CLs
  # use s(..., forcedf=TRUE) to force output to be a data.frame
}
}
