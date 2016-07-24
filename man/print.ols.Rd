\name{print.ols}
\alias{print.ols}
\title{Print ols}
\description{
formatted printing of an object of class \code{ols} using methods taken from
\code{print.lm} and \code{summary.lm}. Prints R-squared, adjusted R-squared,
parameter estimates,
standard errors, and t-statistics (Z statistics if penalized estimation was used).
For penalized
estimation, prints the maximum penalized likelihood estimate of the residual
standard deviation (\code{Sigma}) instead of the usual root mean squared error.
}
\usage{
\method{print}{ols}(x, digits=4, long=FALSE, coefs=TRUE, latex=FALSE, md=FALSE,
title="Linear Regression Model", \dots)
}
\arguments{
\item{x}{fit object}
\item{digits}{number of significant digits to print}
\item{long}{set to \code{TRUE} to print the correlation matrix of
  parameter estimates}
\item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{latex}{a logical value indicating whether information should be
  formatted as plain text or as LaTeX markup}
\item{md}{a logical directing that LaTeX and html be composed for use
	with R Markdown}
\item{title}{a character string title to be passed to \code{prModFit}}
\item{\dots}{other parameters to pass to \code{print} or \code{format}}
}
\seealso{
  \code{\link{ols}}, \code{\link{lm}},\code{\link{prModFit}}
}
\keyword{print}
