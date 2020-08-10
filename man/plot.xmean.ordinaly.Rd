\name{plot.xmean.ordinaly}
\alias{plot.xmean.ordinaly}
\title{
Plot Mean X vs. Ordinal Y
}
\description{
Separately for each predictor variable \eqn{X} in a formula, plots the mean of
\eqn{X} vs. levels of \eqn{Y}.  Then under the proportional odds assumption,
the expected value of the predictor for each \eqn{Y} value is also plotted (as
a dotted line).  This plot is useful for assessing the ordinality assumption 
for \eqn{Y} separately for each \eqn{X}, and for assessing the proportional odds
assumption in a simple univariable way.  If several predictors do not
distinguish adjacent categories of \eqn{Y}, those levels may need to be 
pooled.  This display assumes
that each predictor is linearly related to the log odds of each event in
the proportional odds model.  There is also an option to plot the
expected means assuming a forward continuation ratio model.
}
\usage{
\method{plot}{xmean.ordinaly}(x, data, subset, na.action, subn=TRUE,
                    cr=FALSE, topcats=1, cex.points=.75, \dots)
}
\arguments{
\item{x}{
an S formula.  Response variable is treated as ordinal.  For categorical
predictors, a binary version of the variable is substituted, specifying
whether or not the variable equals the modal category.  Interactions or
non-linear effects are not allowed.
}
\item{data}{
a data frame or frame number
}
\item{subset}{
vector of subscripts or logical vector describing subset of data to
analyze
}
\item{na.action}{
defaults to \code{na.keep} so all NAs are initially retained.  Then NAs
are deleted only for each predictor currently being plotted.
Specify \code{na.action=na.delete} to remove observations that are missing
on any of the predictors (or the response).
}
\item{subn}{
set to \code{FALSE} to suppress a left bottom subtitle specifying the sample size
used in constructing each plot
}
\item{cr}{
set to \code{TRUE} to plot expected values by levels of the response,
assuming a forward continuation ratio model holds.  The function is fairly slow
when this option is specified.
}
\item{topcats}{When a predictor is categorical, by default only the
  proportion of observations in the overall most frequent category will
  be plotted against response variable strata.  Specify a higher value
  of \code{topcats} to make separate plots for the proportion in the
  \code{k} most frequent predictor categories, where \code{k} is
  \code{min(ncat-1, topcats)} and \code{ncat} is the number of unique
  values of the predictor.}
\item{cex.points}{if \code{cr} is \code{TRUE}, specifies the size of the
  \code{"C"} that is plotted.  Default is 0.75.}
\item{...}{
other arguments passed to \code{plot} and \code{lines}
}}
\section{Side Effects}{
plots
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
fh@fharrell.com
}
\references{
Harrell FE et al. (1998): Development of a clinical prediction model for
an ordinal outcome. Stat in Med 17:909--44.
}
\seealso{
\code{\link{lrm}}, \code{\link{residuals.lrm}}, \code{\link{cr.setup}},
\code{\link[Hmisc]{summary.formula}}, \code{\link[Hmisc]{biVar}}.
}
\examples{
# Simulate data from a population proportional odds model
set.seed(1)
n <- 400
age <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
region <- factor(sample(c('north','south','east','west'), n, replace=TRUE))
L <- .2*(age-50) + .1*(blood.pressure-120)
p12 <- plogis(L)    # Pr(Y>=1)
p2  <- plogis(L-1)  # Pr(Y=2)
p   <- cbind(1-p12, p12-p2, p2)   # individual class probabilites
# Cumulative probabilities:
cp  <- matrix(cumsum(t(p)) - rep(0:(n-1), rep(3,n)), byrow=TRUE, ncol=3)
y   <- (cp < runif(n)) \%*\% rep(1,3)
# Thanks to Dave Krantz <dhk@paradox.psych.columbia.edu> for this trick

par(mfrow=c(2,2))
plot.xmean.ordinaly(y ~ age + blood.pressure + region, cr=TRUE, topcats=2)
par(mfrow=c(1,1))
# Note that for unimportant predictors we don't care very much about the
# shapes of these plots.  Use the Hmisc chiSquare function to compute
# Pearson chi-square statistics to rank the variables by unadjusted
# importance without assuming any ordering of the response:
chiSquare(y ~ age + blood.pressure + region, g=3)
chiSquare(y ~ age + blood.pressure + region, g=5)
}
\keyword{category}
\keyword{models}
\keyword{regression}
\keyword{hplot}
\concept{model validation}
\concept{logistic regression model}
