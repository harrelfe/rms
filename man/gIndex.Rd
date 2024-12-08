\name{gIndex}
\alias{gIndex}
\alias{print.gIndex}
\alias{plot.gIndex}
\title{Calculate Total and Partial g-indexes for an rms Fit}
\description{
  \code{gIndex} computes the total \eqn{g}-index for a model based on
  the vector of linear predictors, and the partial \eqn{g}-index for
  each predictor in a model.  The latter is computed by summing all the
  terms involving each variable, weighted by their regression
  coefficients, then computing Gini's mean difference on this sum.  For
  example, a regression model having age and sex and age*sex on the
  right hand side, with corresponding regression coefficients \eqn{b_{1},
	b_{2}, b_{3}}{b1, b2, b3} will have the \eqn{g}-index for age
  computed from Gini's mean
  difference on the product of age \eqn{\times (b_{1} + b_{3}w)}{times
	(b1 + b3*w)} where
  \eqn{w} is an indicator set to one for observations with sex not equal
  to the reference value.  When there are nonlinear terms associated
  with a predictor, these terms will also be combined.

  A \code{print}
  method is defined, and there is a \code{plot} method for displaying
  \eqn{g}-indexes using a dot chart.

	These functions use \code{Hmisc::GiniMd}.
}
\usage{
gIndex(object, partials=TRUE, type=c('ccterms', 'cterms', 'terms'),
           lplabel=if(length(object$scale) && is.character(object$scale))
           object$scale[1] else 'X*Beta',
           fun, funlabel=if(missing(fun)) character(0) else
           deparse(substitute(fun)),
           postfun=if(length(object$scale)==2) exp else NULL,
           postlabel=if(length(postfun))
           ifelse(missing(postfun),
                  if((length(object$scale) > 1) &&
                     is.character(object$scale)) object$scale[2] else
                     'Anti-log',
                     deparse(substitute(postfun))) else character(0),
           \dots)

\method{print}{gIndex}(x, digits=4, abbrev=FALSE,
 vnames=c("names","labels"), \dots)

\method{plot}{gIndex}(x, what=c('pre', 'post'),
 xlab=NULL, pch=16, rm.totals=FALSE,
sort=c('descending', 'ascending', 'none'), \dots)
}
\arguments{
  \item{object}{result of an \code{rms} fitting function}
  \item{partials}{set to \code{FALSE} to suppress computation of partial
	\eqn{g}s}
  \item{type}{defaults to \code{'ccterms'} which causes partial discrimination
	indexes to be computed after maximally combining all related main
	effects and interactions.  The is usually the only way that makes
	sense when considering partial linear predictors.  Specify
	\code{type='cterms'} to only combine a main effect
	with interactions containing it, not also with other main effects
	connected through interactions.  Use \code{type='terms'} to separate
  interactions into their own effects.}
  \item{lplabel}{a replacement for default values such as
	\code{"X*Beta"} or \code{"log odds"}/}
  \item{fun}{an optional function to transform the linear predictors
	before computing the total (only) \eqn{g}.  When this is present, a
	new component \code{gtrans} is added to the attributes of the object
	resulting from \code{gIndex}.}
  \item{funlabel}{a character string label for \code{fun}, otherwise
	taken from the function name itself}
  \item{postfun}{a function to transform \eqn{g} such as \code{exp}
	(anti-log), which is the default for certain models such as the
	logistic and Cox models}
  \item{postlabel}{a label for \code{postfun}}
\item{\dots}{
  For \code{gIndex}, passed to \code{predict.rms}.
  Ignored for \code{print}.  Passed to \code{\link[Hmisc]{dotchart2}}
  for \code{plot}.
}
\item{x}{
  an object created by \code{gIndex} (for \code{print} or \code{plot})
}
\item{digits}{causes rounding to the \code{digits} decimal place}
\item{abbrev}{set to \code{TRUE} to abbreviate labels if
  \code{vname="labels"}}
\item{vnames}{set to \code{"labels"} to print predictor labels instead
  of names}
\item{what}{set to \code{"post"} to plot the transformed \eqn{g}-index
  if there is one (e.g., ratio scale)}
\item{xlab}{\eqn{x}-axis label; constructed by default}
\item{pch}{plotting character for point}
\item{rm.totals}{set to \code{TRUE} to remove the total \eqn{g}-index
  when plotting}
\item{sort}{specifies how to sort predictors by \eqn{g}-index; default
  is in descending order going down the dot chart}
}
\details{
  For stratification factors in a Cox proportional hazards model, there is
  no contribution of variation towards computing a partial \eqn{g}
  except from terms that interact with the stratification variable.
}
\value{
  \code{gIndex} returns a matrix of class \code{"gIndex"} with auxiliary
  information stored as attributes, such as variable labels.
  \code{GiniMd} returns a scalar.
}
\references{
David HA (1968): Gini's mean difference rediscovered.  Biometrika 55:573--575.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
\email{fh@fharrell.com}
}
\seealso{\code{\link{predict.rms}},\code{\link[Hmisc]{GiniMd}}}
\examples{
set.seed(1)
n <- 40
x <- 1:n
w <- factor(sample(c('a','b'), n, TRUE))
u <- factor(sample(c('A','B'), n, TRUE))
y <- .01*x + .2*(w=='b') + .3*(u=='B') + .2*(w=='b' & u=='B') + rnorm(n)/5
dd <- datadist(x,w,u); options(datadist='dd')
f <- ols(y ~ x*w*u, x=TRUE, y=TRUE)
f
anova(f)
z <- list()
for(type in c('terms','cterms','ccterms'))
  {
    zc <- predict(f, type=type)
    cat('type:', type, '\n')
    print(zc)
    z[[type]] <- zc
  }

zc <- z$cterms
GiniMd(zc[, 1])
GiniMd(zc[, 2])
GiniMd(zc[, 3])
GiniMd(f$linear.predictors)
g <- gIndex(f)
g
g['Total',]
gIndex(f, partials=FALSE)
gIndex(f, type='cterms')
gIndex(f, type='terms')

y <- y > .8
f <- lrm(y ~ x * w * u, x=TRUE, y=TRUE, reltol=1e-5)

gIndex(f, fun=plogis, funlabel='Prob[y=1]')

# Manual calculation of combined main effect + interaction effort of
# sex in a 2x2 design with treatments A B, sexes F M,
# model -.1 + .3*(treat=='B') + .5*(sex=='M') + .4*(treat=='B' & sex=='M')

set.seed(1)
X <- expand.grid(treat=c('A','B'), sex=c('F', 'M'))
a <- 3; b <- 7; c <- 13; d <- 5
X <- rbind(X[rep(1, a),], X[rep(2, b),], X[rep(3, c),], X[rep(4, d),])
y <- with(X, -.1 + .3*(treat=='B') + .5*(sex=='M') + .4*(treat=='B' & sex=='M'))
f <- ols(y ~ treat*sex, data=X, x=TRUE)
gIndex(f, type='cterms')
k <- coef(f)
b1 <- k[2]; b2 <- k[3]; b3 <- k[4]
n <- nrow(X)
( (a+b)*c*abs(b2) + (a+b)*d*abs(b2+b3) + c*d*abs(b3))/(n*(n-1)/2 )

# Manual calculation for combined age effect in a model with sex,
# age, and age*sex interaction

a <- 13; b <- 7
sex <- c(rep('female',a), rep('male',b))
agef <- round(runif(a, 20, 30))
agem <- round(runif(b, 20, 40))
age  <- c(agef, agem)
y <- (sex=='male') + age/10 - (sex=='male')*age/20
f <- ols(y ~ sex*age, x=TRUE)
f
gIndex(f, type='cterms')
k <- coef(f)
b1 <- k[2]; b2 <- k[3]; b3 <- k[4]
n <- a + b
sp <- function(w, z=w) sum(outer(w, z, function(u, v) abs(u-v)))

( abs(b2)*sp(agef) + abs(b2+b3)*sp(agem) + 2*sp(b2*agef, (b2+b3)*agem) ) / (n*(n-1))

( abs(b2)*GiniMd(agef)*a*(a-1) + abs(b2+b3)*GiniMd(agem)*b*(b-1) +
  2*sp(b2*agef, (b2+b3)*agem) ) / (n*(n-1))

\dontrun{
# Compare partial and total g-indexes over many random fits
plot(NA, NA, xlim=c(0,3), ylim=c(0,3), xlab='Global',
     ylab='x1 (black)  x2 (red)  x3 (green)  x4 (blue)')
abline(a=0, b=1, col=gray(.9))
big <- integer(3)
n <- 50   # try with n=7 - see lots of exceptions esp. for interacting var
for(i in 1:100) {
   x1 <- runif(n)
   x2 <- runif(n)
   x3 <- runif(n)
   x4 <- runif(n)
   y  <- x1 + x2 + x3 + x4 + 2*runif(n)
   f <- ols(y ~ x1*x2+x3+x4, x=TRUE)
   # f <- ols(y ~ x1+x2+x3+x4, x=TRUE)   # also try this
   w <- gIndex(f)[,1]
   gt <- w['Total']
   points(gt, w['x1, x2'])
   points(gt, w['x3'], col='green')
   points(gt, w['x4'], col='blue')
   big[1] <- big[1] + (w['x1, x2'] > gt)
   big[2] <- big[2] + (w['x3'] > gt)
   big[3] <- big[3] + (w['x4'] > gt)
   }
print(big)
}

options(datadist=NULL)
}
\keyword{predictive accuracy}
\keyword{robust}
\keyword{univar}

