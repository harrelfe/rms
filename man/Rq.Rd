\name{Rq}
\Rdversion{1.1}
\alias{Rq}
\alias{RqFit}
\alias{print.Rq}
\alias{latex.Rq}
\alias{predict.Rq}
\title{rms Package Interface to quantreg Package}
\description{
The \code{Rq} function is the \code{rms} front-end to the
\code{quantreg} package's \code{rq} function.  \code{print} and
\code{latex} methods are also provided, and a fitting function
\code{RqFit} is defined for use in bootstrapping, etc.  Its result is a
function definition.

For the \code{print} method, format of output is controlled by the
user previously running \code{options(prType="lang")} where
\code{lang} is \code{"plain"} (the default), \code{"latex"}, or
\code{"html"}.  For the \code{latex} method, \code{html} will actually
be used of \code{options(prType='html')}.  When using html with Quarto
or RMarkdown, \code{results='asis'} need not be written in the chunk header.
}
\usage{
Rq(formula, tau = 0.5, data=environment(formula),
   subset, weights, na.action=na.delete,
   method = "br", model = FALSE, contrasts = NULL,
   se = "nid", hs = TRUE, x = FALSE, y = FALSE, ...)

\method{print}{Rq}(x, digits=4, coefs=TRUE, title, \dots)

\method{latex}{Rq}(object,
           file = '', append=FALSE,
           which, varnames, columns=65, inline=FALSE, caption=NULL, ...)

\method{predict}{Rq}(object, \dots, kint=1, se.fit=FALSE)

RqFit(fit, wallow=TRUE, passdots=FALSE)
}
\arguments{
  \item{formula}{model formula}
  \item{tau}{
	the single quantile to estimate.  Unlike \code{rq} you cannot estimate
	more than one quantile at one model fitting.
  }
  \item{data,subset,weights,na.action,method,model,contrasts,se,hs}{see
	\code{\link[quantreg]{rq}}}
  \item{x}{set to \code{TRUE} to store the design matrix with the fit.
	For \code{print} is an \code{Rq} object.}
  \item{y}{set to \code{TRUE} to store the response vector with the fit}
  \item{\dots}{
	other arguments passed to one of the \code{rq} fitting routines.
	For \code{latex.Rq} these are optional arguments passed to
	\code{latexrms}.  Ignored for \code{print.Rq}.  For
	\code{predict.Rq} this is usually just a \code{newdata} argument.
  }
  \item{digits}{
	number of significant digits used in formatting results in
	\code{print.Rq}.
  }
  \item{coefs}{specify \code{coefs=FALSE} to suppress printing the table
	of model coefficients, standard errors, etc.  Specify \code{coefs=n}
	to print only the first \code{n} regression coefficients in the
	model.}
  \item{title}{a character string title to be passed to \code{prModFit}}
  \item{object}{an object created by \code{Rq}}
  \item{file,append,which,varnames,columns,inline,caption}{see
		\code{\link{latexrms}}}
	\item{kint}{ignored}
  \item{se.fit}{set to \code{TRUE} to obtain standard errors of
	predicted quantiles}
  \item{fit}{an object created by \code{Rq}}
  \item{wallow}{
	set to \code{TRUE} if \code{weights} are allowed in the
	current context.
  }
  \item{passdots}{
	set to \code{TRUE} if \dots may be passed to the fitter}
}
\value{
  \code{Rq} returns a list of class \code{"rms", "lassorq"} or \code{"scadrq",
	"Rq"}, and \code{"rq"}.  \code{RqFit} returns a function
  definition.  \code{latex.Rq} returns an object of class \code{"latex"}.
}
\author{
  Frank Harrell
}
\note{
  The author and developer of methodology in the \code{quantreg} package
  is Roger Koenker.
}
\seealso{
  \code{\link[quantreg]{rq}}, \code{\link{prModFit}}, \code{\link{orm}}
}
\examples{
\dontrun{
set.seed(1)
n <- 100
x1 <- rnorm(n)
y <- exp(x1 + rnorm(n)/4)
dd <- datadist(x1); options(datadist='dd')
fq2 <- Rq(y ~ pol(x1,2))
anova(fq2)
fq3 <- Rq(y ~ pol(x1,2), tau=.75)
anova(fq3)
pq2 <- Predict(fq2, x1)
pq3 <- Predict(fq3, x1)
p <- rbind(Median=pq2, Q3=pq3)
plot(p, ~ x1 | .set.)
# For superpositioning, with true curves superimposed
a <- function(x, y, ...) {
 x <- unique(x)
 col <- trellis.par.get('superpose.line')$col
 llines(x, exp(x), col=col[1], lty=2)
 llines(x, exp(x + qnorm(.75)/4), col=col[2], lty=2)
}
plot(p, addpanel=a)
}
}
\keyword{models}
\keyword{nonparametric}
