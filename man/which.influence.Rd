\name{which.influence}
\alias{which.influence}
\alias{show.influence}
\title{
Which Observations are Influential
}
\description{
Creates a list with a component for each factor in the model.  The
names of the components are the factor names.  Each component contains
the observation identifiers of all observations that are "overly
influential" with respect to that factor, meaning that \eqn{|dfbetas| >
u} for at least one \eqn{\beta_i}{beta i} associated with that factor,
for a given \code{cutoff}.  The default \code{cutoff} is \code{.2}.  The
fit must come from a function that has \code{resid(fit, type="dfbetas")}
defined.

\code{show.influence}, written by Jens Oehlschlaegel-Akiyoshi, applies the
result of \code{which.influence} to a data frame, usually the one used to
fit the model, to report the results.
}
\usage{
which.influence(fit, cutoff=.2)

show.influence(object, dframe, report=NULL, sig=NULL, id=NULL)
}
\arguments{
\item{fit}{
fit object
}
\item{object}{
the result of \code{which.influence}
}
\item{dframe}{
data frame containing observations pertinent to the model fit
}
\item{cutoff}{
cutoff value
}
\item{report}{
other columns of the data frame to report besides those corresponding
to predictors that are influential for some observations
}
\item{sig}{
runs results through \code{signif} with \code{sig} digits if \code{sig} is given
}
\item{id}{
a character vector that labels rows of \code{dframe} if \code{row.names} were
not used
}}
\value{
\code{show.influence} returns a marked dataframe with the first column being
a count of influence values
}
\author{
Frank Harrell\cr
Department of Biostatistics, Vanderbilt University\cr
fh@fharrell.com
\cr

Jens Oehlschlaegel-Akiyoshi\cr
Center for Psychotherapy Research\cr
Christian-Belser-Strasse 79a\cr
D-70597 Stuttgart Germany\cr
oehl@psyres-stuttgart.de
}
\seealso{
\code{\link{residuals.lrm}}, \code{\link{residuals.cph}},
\code{\link{residuals.ols}}, \code{\link{rms}}, \code{\link{lrm}},
\code{\link{ols}}, \code{\link{cph}}
}
\examples{
#print observations in data frame that are influential,
#separately for each factor in the model
x1 <- 1:20
x2 <- abs(x1-10)
x3 <- factor(rep(0:2,length.out=20))
y  <- c(rep(0:1,8),1,1,1,1)
f  <- lrm(y ~ rcs(x1,3) + x2 + x3, x=TRUE,y=TRUE)
w <- which.influence(f, .55)
nam <- names(w)
d   <- data.frame(x1,x2,x3,y)
for(i in 1:length(nam)) {
 print(paste("Influential observations for effect of ",nam[i]),quote=FALSE)
 print(d[w[[i]],])
}

show.influence(w, d)  # better way to show results
}
\keyword{models}
\keyword{regression}
\keyword{survival}
\concept{logistic regression model}
