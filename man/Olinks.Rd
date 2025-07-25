% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Olinks.r
\name{Olinks}
\alias{Olinks}
\title{Likehood-Based Statistics for Other Links for orm Fits}
\usage{
Olinks(
  object,
  links = c("logistic", "probit", "loglog", "cloglog"),
  dec = 3,
  gradtol = 0.001
)
}
\arguments{
\item{object}{an object created by \code{orm} with \verb{x=TRUE, y=TRUE}}

\item{links}{a vector of links to consider other than the one used to get \code{object}}

\item{dec}{number of digits to the right of the decimal place to round statistics to}

\item{...}{passed to \code{\link[=orm.fit]{orm.fit()}}}
}
\value{
data frame.  The \code{R2} column is from the last adjusted \eqn{R^2} computed by \code{orm},
which adjustes for the effective sample size and the number of betas.
}
\description{
Likehood-Based Statistics for Other Links for orm Fits
}
\examples{
\dontrun{
f <- orm(y ~ x1 + x2, family='loglog', x=TRUE, y=TRUE)
Olinks(f)
}
}
\author{
Frank Harrell
}
