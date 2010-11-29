\name{specs.rms}
\alias{specs.rms}
\alias{specs}
\alias{print.specs.rms}
\title{rms Specifications for Models}
\description{
Prints the design specifications, e.g., number of parameters for each
factor, levels of categorical factors, knot locations in splines,
pre-transformations, etc. 
}
\usage{
specs(fit, \dots)
\method{specs}{rms}(fit, long=FALSE, \dots)

\method{print}{specs.rms}(x, \dots)
}
\arguments{
\item{fit}{
a fit object created with the \code{rms} library in effect
}
\item{x}{
an object returned by \code{specs}
}
\item{long}{
if  \code{TRUE}, causes the plotting and estimation limits to be printed for
each factor
}
\item{\dots}{ignored}
}
\value{
a list containing information about the fit and the predictors as elements
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\seealso{
\code{\link{rms}}, \code{\link{rms.trans}}, \code{\link{latex.rms}}, \code{\link{datadist}}
}
\examples{
set.seed(1)
blood.pressure <- rnorm(200, 120, 15)
dd <- datadist(blood.pressure)
options(datadist='dd')
L <- .03*(blood.pressure-120)
sick <- ifelse(runif(200) <= plogis(L), 1, 0)
f <- lrm(sick ~ rcs(blood.pressure,5))
specs(f)    # find out where 5 knots are placed
g <- Glm(sick ~ rcs(blood.pressure,5), family=binomial)
specs(g,long=TRUE)
options(datadist=NULL)
}
\keyword{models}
\keyword{regression}
\keyword{methods}
