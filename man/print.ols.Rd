\name{print.ols}
\alias{print.ols}
\title{Print ols}
\description{
Formatted printing of an object of class \code{ols} using methods taken
from \code{print.lm} and \code{summary.lm}. Prints R-squared, adjusted
R-squared, parameter estimates, standard errors, and t-statistics (Z
statistics if penalized estimation was used).  For penalized estimation,
prints the maximum penalized likelihood estimate of the residual
standard deviation (\code{Sigma}) instead of the usual root mean squared
error.
Format of output is controlled by the user previously running
\code{options(prType="lang")} where \code{lang} is \code{"plain"} (the default),
\code{"latex"}, or \code{"html"}.  When using html with Quarto or RMarkdown,
\code{results='asis'} need not be written in the chunk header.
}
\usage{
\method{print}{ols}(x, digits=4, long=FALSE, coefs=TRUE, 
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
\item{title}{a character string title to be passed to \code{prModFit}}
\item{\dots}{other parameters to pass to \code{print} or \code{format}}
}
\seealso{
  \code{\link{ols}}, \code{\link{lm}},\code{\link{prModFit}}
}
\keyword{print}
