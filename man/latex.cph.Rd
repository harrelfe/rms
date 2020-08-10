\name{latex.cph}
\alias{latex.cph}
\alias{latex.lrm}
\alias{latex.ols}
\alias{latex.orm}
\alias{latex.pphsm}
\alias{latex.psm}
\title{LaTeX Representation of a Fitted Cox Model}

\description{Creates a file containing a LaTeX representation of the fitted model.}

\usage{
\method{latex}{cph}(object, title,
      file='', 
      append=FALSE, surv=TRUE, maxt=FALSE, which=NULL, varnames, columns=65, 
      inline=FALSE, before=if(inline)"" else "& &", after="", dec=3,
      pretrans=TRUE, caption, digits=.Options$digits, size="",
      \dots) # for cph fit

\method{latex}{lrm}(object, title, file, append, which, varnames,
columns, inline, before, after, pretrans, caption,
digits=.Options$digits, size="", \dots) # for lrm fit

\method{latex}{ols}(object, title, file, append, which, varnames,
columns, inline, before, after, pretrans, caption,
digits=.Options$digits, size="", \dots) # ols fit

\method{latex}{orm}(object, title, file, append, which, varnames,
columns, inline, before, after, pretrans, caption,
digits=.Options$digits, size="", intercepts=nrp < 10, \dots) # for orm fit

\method{latex}{pphsm}(object, title, file, append, which=NULL, varnames,
columns, inline, before, after, pretrans, caption,
digits=.Options$digits, size="", \dots) # pphsm fit

\method{latex}{psm}(object, title, file, append, which=NULL, varnames,
columns, inline, before, after, pretrans, caption,
digits=.Options$digits, size="", \dots) # psm fit
}
\arguments{
\item{object}{
a fit object created by a \code{rms} fitting function.
}
\item{title}{ignored}
\item{file,append}{see \code{\link[Hmisc:latex]{latex.default}}.  Defaults to
	the console.  When usinghtml/markdown, \code{file} is ignored.}
\item{surv}{
if \code{surv=TRUE} was specified to \code{cph}, the underlying survival
probabilities from \code{object$surv.summary} will be placed in a table
unless \code{surv=FALSE}.
}
\item{maxt}{
if the maximum follow-up time in the data (\code{object$maxtime}) exceeds the
last entry in \code{object$surv.summary}, underlying survival estimates at
\code{object$maxtime} will be added to the table if \code{maxt=TRUE}.
}
\item{which,varnames,columns,inline,before,dec,pretrans}{see
  \code{\link[Hmisc]{latex.default}}}
\item{after}{if not an empty string, added to end of markup if
  \code{inline=TRUE}}
\item{caption}{a character string specifying a title for the equation to
  be centered and typeset in bold face.   Default is no title.
}
\item{digits}{see \link{latexrms}}
\item{size}{a LaTeX size to use, without the slash.  Default is the
  prevailing size}
\item{intercepts}{for \code{orm} fits.  Default is to print intercepts
	if they are fewer than 10 in number.  Set to \code{TRUE} or
	\code{FALSE} to force.}
\item{\dots}{ignored}
}
\value{
the name of the created file, with class \code{c("latex","file")}.  This
object works with latex viewing and printing commands in Hmisc.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
}
\seealso{
\code{\link{latexrms}}, \code{\link[Hmisc]{rcspline.restate}},
\code{\link[Hmisc]{latex}}
}
\examples{
\dontrun{
units(ftime) <- "Day"
f <- cph(Surv(ftime, death) ~ rcs(age)+sex, surv=TRUE, time.inc=60)
w <- latex(f, file='f.tex')  #Interprets fitted model and makes table of S0(t)
               #for t=0,60,120,180,...
w              #displays image, if viewer installed and file given above
latex(f)   # send LaTeX code to the console for knitr
options(prType='html')
latex(f)       # for use with knitr and R Markdown using MathJax
}
}
\keyword{regression}
\keyword{character}
\keyword{survival}
\keyword{interface}
\keyword{models}
