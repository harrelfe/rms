\name{ExProb}
\alias{ExProb}
\alias{ExProb.orm}
\alias{plot.ExProb}
\title{Function Generator For Exceedance Probabilities}
\description{
	For an \code{orm} object generates a function for computing the
	estimates of the function Prob(Y>=y) given one or more values of the
	linear predictor using the reference (median) intercept.  This
	function can optionally be evaluated at only a set of user-specified
	\code{y} values, otherwise a right-step function is returned.  There
	is a plot method for plotting the step functions, and if more than one
	linear predictor was evaluated multiple step functions are drawn.
	\code{ExProb} is especially useful for \code{\link{nomogram}}.
	
	Optionally a normal approximation for a confidence
	interval for exceedance probabilities will be computed using the delta
	method, if 
	\code{conf.int > 0} is specified to the function generated from calling
	\code{ExProb}.  In that case, a \code{"lims"} attribute is included
	in the result computed by the derived cumulative probability function.
}
\usage{
ExProb(object, \dots)

\method{ExProb}{orm}(object, codes = FALSE, ...)

\method{plot}{ExProb}(x, \dots, data=NULL,
                      xlim=NULL, xlab=x$yname, ylab=expression(Prob(Y>=y)),
                      col=par('col'), col.vert='gray85', pch=20,
                      pch.data=21, lwd=par('lwd'), lwd.data=lwd,
                      lty.data=2, key=TRUE)
}
\arguments{
  \item{object}{a fit object from \code{orm}}
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
	\code{ExProb} returns an R function.  Running the function returns an
  object of class \code{"ExProb"}.
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
}
