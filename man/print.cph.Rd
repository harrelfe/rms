\name{print.cph}
\alias{print.cph}
\title{Print cph Results}

\description{
Formatted printing of an object of class \code{cph}. Prints strata
frequencies, parameter estimates, standard errors, z-statistics, numbers
of missing values, etc.
Format of output is controlled by the user previously running
\code{options(prType="lang")} where \code{lang} is \code{"plain"} (the default),
\code{"latex"}, or \code{"html"}.  This does not require \code{results='asis'}
in \code{knitr} chunk headers.
}
\usage{
\method{print}{cph}(x, digits=4, r2=c(0,2,4), table=TRUE, conf.int=FALSE, 
coefs=TRUE, pg=FALSE, title='Cox Proportional Hazards Model', \dots)
}
\arguments{
\item{x}{fit object}
\item{digits}{number of digits to right of decimal place to print}
\item{r2}{vector of integers specifying which R^2 measures to print,
		with 0 for Nagelkerke R^2 and 1:4 corresponding to the 4 measures
		computed by \code{\link[Hmisc]{R2Measures}}.  Default is to print
		Nagelkerke (labeled R2) and second and fourth \code{R2Measures}
		which are the measures adjusted for the number of predictors, first
		for the raw sample size then for the effective sample size, which
		here is the number of non-censored observations.}
\item{conf.int}{
set to e.g. .95 to print 0.95 confidence intervals on simple hazard
ratios (which are usually meaningless as one-unit changes are seldom
relevant and most models contain multiple terms per predictor)
}
\item{table}{
set to \code{FALSE} to suppress event frequency statistics
}
\item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
  of model coefficients, standard errors, etc.  Specify \code{coefs=n}
  to print only the first \code{n} regression coefficients in the
  model.}
\item{pg}{set to \code{TRUE} to print g-indexes}
\item{title}{a character string title to be passed to \code{prModFit}}
\item{\dots}{arguments passed to \code{prModFit}}
}
\seealso{
\code{\link[survival]{coxph}}, \code{\link{prModFit}}
}
\keyword{print}
