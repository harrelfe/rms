\name{Glm}
\alias{Glm}
\alias{print.Glm}
\title{rms Version of glm}
\description{
This function saves \code{rms} attributes with the fit object so that
\code{anova.rms}, \code{Predict}, etc. can be used just as with
\code{ols} and other fits.  No \code{validate} or \code{calibrate}
methods exist for \code{Glm} though.
}
\usage{
Glm(formula, family = gaussian, data = list(), weights = NULL, subset =
NULL, na.action = na.delete, start = NULL, offset = NULL, control =
glm.control(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE,
contrasts = NULL, \dots)

\method{print}{Glm}(x, digits=4, coefs=TRUE, latex=FALSE,
title='General Linear Model', \dots)
}
\arguments{
  \item{formula,family,data,weights,subset,na.action,start,offset,control,model,method,x,y,contrasts}{see \code{\link{glm}}; for \code{print}, \code{x} is
	the result of \code{Glm}}
  \item{\dots}{ignored}
  \item{digits}{number of significant digits to print}
  \item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
	of model coefficients, standard errors, etc.  Specify \code{coefs=n}
	to print only the first \code{n} regression coefficients in the
	model.}
  \item{latex}{a logical value indicating whether information should be
	formatted as plain text or as LaTeX markup}
  \item{title}{a character string title to be passed to \code{prModFit}}
}
\value{a fit object like that produced by \code{\link{glm}} but with
  \code{rms} attributes and a \code{class} of \code{"rms"},
  \code{"Glm"}, \code{"glm"}, and \code{"lm"}.  The \code{g}
  element of the fit object is the \eqn{g}-index.
}
\seealso{\code{\link{glm}},\code{\link{rms}},\code{\link{GiniMd}},
   \code{\link{prModFit}},\code{\link[stats]{residuals.glm}}}
\examples{
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
f <- glm(counts ~ outcome + treatment, family=poisson())
f
anova(f)
summary(f)
f <- Glm(counts ~ outcome + treatment, family=poisson())
# could have had rcs( ) etc. if there were continuous predictors
f
anova(f)
summary(f, outcome=c('1','2','3'), treatment=c('1','2','3'))
}
\keyword{models}
\keyword{regression}
