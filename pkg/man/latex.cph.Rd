\name{latex.cph}
\alias{latex.cph}
\alias{latex.lrm}
\alias{latex.ols}
\alias{latex.pphsm}
\alias{latex.psm}
\title{LaTeX Representation of a Fitted Cox Model}

\description{
  Creates a file containing a LaTeX representation of the fitted model. 
}

\usage{
\method{latex}{cph}(object, title,
      file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
      append=FALSE, surv=TRUE, maxt=FALSE, which=NULL, varnames, columns=65, 
      inline=FALSE, before=if(inline)"" else "& &", dec=3,
      pretrans=TRUE, caption, \dots) # for cph fit

\method{latex}{lrm}(object, title, file, append, which, varnames,
columns, inline, before, pretrans, caption, \dots) # for lrm fit

\method{latex}{ols}(object, title, file, append, which, varnames,
columns, inline, before, pretrans, caption, \dots) # ols fit

\method{latex}{pphsm}(object, title, file, append, which=NULL, varnames,
columns, inline, before, pretrans, caption, \dots) # pphsm fit

\method{latex}{psm}(object, title, file, append, which=NULL, varnames,
columns, inline, before, pretrans, caption, \dots) # psm fit
}
\arguments{
\item{object}{
a fit object created by a \code{rms} fitting function.
}
\item{title}{ignored}
\item{file,append}{see \code{\link[Hmisc]{latex.default}}}
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
\item{which,varnames,columns,inline,before,dec,pretrans}{see \code{\link[Hmisc]{latex.default}}}
\item{caption}{a character string specifying a title for the equation to
  be centered and typeset in bold face.   Default is no title.
  }
\item{\dots}{ignored}
}
\value{
the name of the created file, with class \code{c("latex","file")}.  This
object works with latex viewing and printing commands in Hmisc.
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{latex.rms}}, \code{\link[Hmisc]{rcspline.restate}}, \code{\link[Hmisc]{latex}}
}
\examples{
\dontrun{
units(ftime) <- "Day"
f <- cph(Surv(ftime, death) ~ rcs(age)+sex, surv=TRUE, time.inc=60)
w <- latex(f)  #Interprets fitted model and makes table of S0(t)
               #for t=0,60,120,180,\dots   Creates file f.tex
w              #displays image, if viewer installed
latex(f,file="")   # send LaTeX code to the screen
}
}
\keyword{regression}
\keyword{character}
\keyword{survival}
\keyword{interface}
\keyword{models}
