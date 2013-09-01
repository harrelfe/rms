\name{latexrms}
\alias{latexrms}
\alias{latex.bj}
\alias{latex.Glm}
\alias{latex.Gls}
\title{LaTeX Representation of a Fitted Model}
\description{
Creates a file containing a LaTeX representation of the fitted model.  For
model-specific typesetting there is \code{latex.lrm}, \code{latex.cph},
\code{latex.psm} and \code{latex.ols}. \code{latex.cph} has some
arguments that are specific to \code{cph} models.  These routines work
with the \code{display} package from statlib to display and print the
formatted model fits.  \code{latexrms} is the core function which is
called internally by \code{latexrms} (which is called by
\code{latex.cph}, \code{latex.ols}, etc.).
}
\usage{
latexrms(object,
 file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
 append=FALSE, which=1:p, varnames, columns=65, prefix=NULL, inline=FALSE,
 before=if(inline)"" else "& &", after="", intercept, pretrans=TRUE,
 digits=.Options$digits, size="")
}
\arguments{
\item{object}{
a fit object created by a fitting function in the \code{rms} series
}
\item{file}{
  name of \code{.tex} file to create, default is first word of argument
  \code{object} with \code{".tex"} added.  Set to \code{""} to send LaTeX
  output to standard output.
}
\item{append}{whether or not to append to an existing file}
\item{which}{
  a vector of subcripts (corresponding to \code{object$Design$name})
  specifying a submodel to print. Default is to describe the whole
  model. 
  \code{which} can also be a vector of character strings specifying the
  factor names to print. Enough of each string is needed to ensure
  a unique match. Names for interaction effects are of the form
  \code{"age * sex"}. For any interaction effect for which you do not
  request main effects, the main effects will be added to \code{which}.
  When \code{which} is given, the model structural statement is not
  included. In this case, intercepts are not included either.
}
\item{varnames}{
  variable names to substitute for non-interactions. Order must correspond
  to \code{object$Design$name} and interactions must be omitted.
  Default is
  \code{object$Design$name[object$Design$assume.code!=9]}. \code{varnames}
 can contain any LaTeX commands such as subscripts and "\\\\\\\\frac"    
 (all "\\" must be quadrupled.)
 Any "/" must be preceeded by "\\\\" (2, not 4 backslashes).
 Elements of \code{varnames} for interactions are ignored; they can be
 set to any value.
  }
\item{columns}{
  maximum number of columns of printing characters to allow before
  outputting a LaTeX newline command
}
\item{prefix}{
  if given, a LaTeX \\lefteqn command of the form \code{\\lefteqn\{prefix =\} \\\\}
  will be inserted to print a left-hand-side of the equation.
}
\item{inline}{
  Set to \code{TRUE} to create text for insertion in an in-line equation. This
  text contains only the expansion of X beta, and is not surrounded by
  \code{"$"}.
}
\item{before}{
  a character string to place before each line of output. Use the default
  for a LaTeX \code{eqnarray} environment.  For \code{inline=TRUE}, the
 \code{before} string, if not an empty string, will be placed once
 before the entire markup.
}
\item{after}{
  a character string to place after the output if \code{inline=TRUE}
  }
\item{intercept}{
  a special intercept value to include that is not part of the standard
  model parameters (e.g., centering constant in Cox model). Only allowed
  in the \code{latexrms} rendition.
}
\item{pretrans}{
  if any spline or polynomial-expanded variables are themselves
  transformed, a table of pre-transformations will be formed unless
  \code{pretrans=FALSE}.
}
\item{digits}{number of digits of precision to use in formatting
  coefficients and other numbers}
\item{size}{a LaTeX font size to use for the output, without the slash.
 Default is current size.}
}
\value{a file name of class \code{"latex"}}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link[Hmisc]{latex}}, \code{\link[Hmisc]{rcspline.restate}},
\code{\link{rms}}
}
\examples{
\dontrun{
f <- lrm(death ~ rcs(age)+sex)
w <- latex(f)
w     # displays, using e.g. xdvi
latex(f, file="")    # send LaTeX code to screen
}
}
\keyword{models}
\keyword{regression}
\keyword{character}
\keyword{methods}
\keyword{interface}
