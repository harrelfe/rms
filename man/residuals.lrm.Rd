\name{residuals.lrm}
\alias{residuals.lrm}
\alias{plot.lrm.partial}
\title{Residuals from a Logistic Regression Model Fit}
\description{
For a binary logistic model fit, computes the following residuals, letting
\eqn{P} denote the predicted probability of the higher category of \eqn{Y},
\eqn{X} denote the design matrix (with a column of 1s for the intercept), and
\eqn{L} denote the logit or linear predictors: ordinary or Li-Shepherd
(\eqn{Y-P}), score 
(\eqn{X (Y-P)}), pearson (\eqn{(Y-P)/\sqrt{P(1-P)}}), deviance (for \eqn{Y=0} is
\eqn{-\sqrt{2|\log(1-P)|}}, for \eqn{Y=1} is \eqn{\sqrt{2|\log(P)|}},
pseudo dependent variable used in influence statistics 
(\eqn{L + (Y-P)/(P(1-P))}), and partial (\eqn{X_{i}\beta_{i} +
  (Y-P)/(P(1-P))}). 


Will compute all these residuals for an ordinal logistic model, using
as temporary binary responses dichotomizations of \eqn{Y}, along with the corresponding
\eqn{P}, the probability that \eqn{Y \geq} cutoff.  For
\code{type="partial"}, all 
possible dichotomizations are used, and for \code{type="score"}, the actual
components of the first derivative of the log likelihood are used for
an ordinal model.  For \code{type="li.shepherd"} the residual is
\eqn{Pr(W<Y) - Pr(W>Y)} where Y is the observed response and W is a
random variable from the fitted distribution.
Alternatively, specify \code{type="score.binary"}
to use binary model score residuals but for all cutpoints of \eqn{Y}
(plotted only, not returned). The \code{score.binary}, 
\code{partial}, and perhaps \code{score} residuals are useful for checking the proportional odds assumption.
If the option \code{pl=TRUE} is used to plot the \code{score} or \code{score.binary} residuals, 
a score residual plot is
made for each column of the design (predictor) matrix, with \code{Y} cutoffs on the
x-axis and the mean +- 1.96 standard errors of the score residuals on
the y-axis.  You can instead use a box plot to display these residuals,
for both \code{score.binary} and \code{score}.
Proportional odds dictates a horizontal \code{score.binary} plot.  Partial
residual plots use smooth nonparametric estimates, separately for each
cutoff of \eqn{Y}.  One examines that plot for parallelism of the curves
to check the proportional odds assumption, as
well as to see if the predictor behaves linearly.


Also computes a variety of influence statistics and the 
le Cessie - van Houwelingen - Copas - Hosmer unweighted sum of squares test
for global goodness of fit, done separately for each cutoff of \eqn{Y} in the
case of an ordinal model.


The \code{plot.lrm.partial} function computes partial residuals for a series
of binary logistic model fits that all used the same predictors and that
specified \code{x=TRUE, y=TRUE}.  It then computes smoothed partial residual
relationships (using \code{lowess} with \code{iter=0}) and plots them separately
for each predictor, with residual plots from all model fits shown on the
same plot for that predictor.
}
\usage{
\method{residuals}{lrm}(object, type=c("li.shepherd","ordinary",
 "score", "score.binary", "pearson", "deviance", "pseudo.dep",
 "partial", "dfbeta", "dfbetas", "dffit", "dffits", "hat", "gof", "lp1"),
           pl=FALSE, xlim, ylim, kint, label.curves=TRUE, which, \dots)

\method{plot}{lrm.partial}(\dots, labels, center=FALSE, ylim)
}
\arguments{
\item{object}{
object created by \code{lrm}
}
\item{\dots}{
for \code{residuals}, applies to \code{type="partial"}  when \code{pl}
is not \code{FALSE}.  These are extra arguments passed to the smoothing
function.  Can also be used to pass extra arguments to \code{boxplot}
for \code{type="score"} or \code{"score.binary"}. 
For \code{plot.lrm.partial} this specifies a series of binary model fit
objects.
}
\item{type}{
type of residual desired.  Use \code{type="lp1"} to get approximate leave-out-1
linear predictors, derived by subtracting the \code{dffit} from the original
linear predictor values.
}
\item{pl}{
applies only to \code{type="partial"}, \code{"score"}, and \code{"score.binary"}.  
For score residuals in an ordinal model, set \code{pl=TRUE} to get means and 
approximate 0.95 confidence bars
vs. \eqn{Y}, separately for each \eqn{X}.  Alternatively, specify
\code{pl="boxplot"} to use \code{boxplot} to
draw the plot, with notches and with width proportional to the square
root of the cell sizes.
For partial residuals, set \code{pl=TRUE} (which uses \code{lowess}) or \code{pl="supsmu"}
to get smoothed partial
residual plots for all columns of \eqn{X} using \code{supsmu}.
Use \code{pl="loess"} to use \code{loess} and get confidence bands (\code{"loess"} is not
implemented for ordinal responses).  Under R, \code{pl="loess"} uses
\code{lowess} and does not provide confidence bands.
If there is more than one \eqn{X},
you should probably use \code{par(mfrow=c( , ))} before calling \code{resid}.
Note that \code{pl="loess"} results in \code{plot.loess} being called, which
requires a large memory allocation.
}
\item{xlim}{
plotting range for x-axis (default = whole range of predictor)
}
\item{ylim}{
plotting range for y-axis (default = whole range of residuals,
range of all confidence intervals for \code{score} or \code{score.binary} or range
of all smoothed curves for \code{partial} if
\code{pl=TRUE}, or 0.1 and 0.9 quantiles of the residuals for \code{pl="boxplot"}.)
}
\item{kint}{
for an ordinal model for residuals other than \code{li.shepherd}, \code{partial}, \code{score}, or
\code{score.binary}, specifies
the intercept (and the cutoff of \eqn{Y}) to use for the calculations.
Specifying \code{kint=2}, for example, means to use \eqn{Y \geq} 3rd level.
}
\item{label.curves}{
set to \code{FALSE} to suppress curve labels when \code{type="partial"}.  The default,
\code{TRUE}, causes \code{labcurve} to be invoked to label curves where they are most
separated.  \code{label.curves} can be a list containing the \code{opts} parameter
for \code{labcurve}, to send options to \code{labcurve}, such as \code{tilt}.  The
default for \code{tilt} here is \code{TRUE}.
}
\item{which}{
a vector of integers specifying column numbers of the design matrix for which to compute or plot residuals, for \code{type="partial","score","score.binary"}.
}
\item{labels}{
for \code{plot.lrm.partial} this specifies a vector of character strings 
providing labels for the list of binary fits.  By default, the names of
the fit objects are used as labels.  The \code{labcurve} function is used
to label the curve with the \code{labels}.
}
\item{center}{
for \code{plot.lrm.partial} this causes partial residuals for every model to have a mean of zero before smoothing and plotting
}}
\value{
a matrix (\code{type="partial","dfbeta","dfbetas","score"}), 
test statistic (\code{type="gof"}), or a vector otherwise.  
For partial residuals from an ordinal
model, the returned object is a 3-way array (rows of \eqn{X} by columns
of \eqn{X} by cutoffs of \eqn{Y}), and NAs deleted during the fit
are not re-inserted into the residuals.  For \code{score.binary}, nothing
is returned.
}
\details{
For the goodness-of-fit test, the le Cessie-van Houwelingen normal test
statistic for the unweighted sum of squared errors (Brier score times \eqn{n})
is used.  For an ordinal response variable, the test 
for predicting the probability that \eqn{Y\geq j} is done separately for
all \eqn{j} (except the first).  Note that the test statistic can have strange behavior 
(i.e., it is far too large) if the model has no predictive value.


For most of the values of \code{type}, you must have specified \code{x=TRUE, y=TRUE} to
\code{lrm}.


There is yet no literature on interpreting score residual plots for the
ordinal model.  Simulations when proportional odds is satisfied have
still shown a U-shaped residual plot.  The series of binary model score
residuals for all cutoffs of \eqn{Y} seems to better check the assumptions.
See the examples.

The li.shepherd residual is a single value per observation on the
probability scale and can be useful for examining linearity, checking
for outliers, and measuring residual correlation.
}
\author{
Frank Harrell\cr
Department of Biostatistics\cr
Vanderbilt University\cr
f.harrell@vanderbilt.edu
}
\references{
Landwehr, Pregibon, Shoemaker. JASA 79:61--83, 1984.


le Cessie S, van Houwelingen JC. Biometrics 47:1267--1282, 1991.


Hosmer DW, Hosmer T, Lemeshow S, le Cessie S, Lemeshow S.  A
comparison of goodness-of-fit tests for the logistic regression model.
Stat in Med 16:965--980, 1997.


Copas JB.  Applied Statistics 38:71--80, 1989.

Li C, Shepherd BE.  Biometrika 99:473-480, 2012.
}
\seealso{
\code{\link{lrm}}, \code{\link{naresid}}, \code{\link{which.influence}},
\code{\link{loess}}, \code{\link{supsmu}}, \code{\link{lowess}},
\code{\link{boxplot}}, \code{\link[Hmisc]{labcurve}}
}
\examples{
set.seed(1)
x1 <- runif(200, -1, 1)
x2 <- runif(200, -1, 1)
L  <- x1^2 - .5 + x2
y  <- ifelse(runif(200) <= plogis(L), 1, 0)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
resid(f)            #add rows for NAs back to data
resid(f, "score")   #also adds back rows
r <- resid(f, "partial")  #for checking transformations of X's
par(mfrow=c(1,2))
for(i in 1:2) {
  xx <- if(i==1)x1 else x2
  plot(xx, r[,i], xlab=c('x1','x2')[i])
  lines(lowess(xx,r[,i]))
}
resid(f, "partial", pl="loess")  #same as last 3 lines
resid(f, "partial", pl=TRUE) #plots for all columns of X using supsmu
resid(f, "gof")           #global test of goodness of fit
lp1 <- resid(f, "lp1")    #approx. leave-out-1 linear predictors
-2*sum(y*lp1 + log(1-plogis(lp1)))  #approx leave-out-1 deviance
                                    #formula assumes y is binary


# Simulate data from a population proportional odds model
set.seed(1)
n   <- 400
age <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
L <- .05*(age-50) + .03*(blood.pressure-120)
p12 <- plogis(L)    # Pr(Y>=1)
p2  <- plogis(L-1)  # Pr(Y=2)
p   <- cbind(1-p12, p12-p2, p2)   # individual class probabilites
# Cumulative probabilities:
cp  <- matrix(cumsum(t(p)) - rep(0:(n-1), rep(3,n)), byrow=TRUE, ncol=3)
# simulate multinomial with varying probs:
y <- (cp < runif(n)) \%*\% rep(1,3)
y <- as.vector(y)
# Thanks to Dave Krantz for this trick
f <- lrm(y ~ age + blood.pressure, x=TRUE, y=TRUE)
par(mfrow=c(2,2))
resid(f, 'score.binary',   pl=TRUE)              #plot score residuals
resid(f, 'partial', pl=TRUE)                     #plot partial residuals
resid(f, 'gof')           #test GOF for each level separately


# Show use of Li-Shepherd residuals
f.wrong <- lrm(y ~ blood.pressure, x=TRUE, y=TRUE)
par(mfrow=c(2,1))
# li.shepherd residuals from model without age
plot(age, resid(f.wrong, type="li.shepherd"),
     ylab="li.shepherd residual")
lines(lowess(age, resid(f.wrong, type="li.shepherd")))
# li.shepherd residuals from model including age
plot(age, resid(f, type="li.shepherd"),
     ylab="li.shepherd residual")
lines(lowess(age, resid(f, type="li.shepherd")))


# Make a series of binary fits and draw 2 partial residual plots
#
f1 <- lrm(y>=1 ~ age + blood.pressure, x=TRUE, y=TRUE)
f2  <- update(f1, y==2 ~.)
par(mfrow=c(2,1))
plot.lrm.partial(f1, f2)


# Simulate data from both a proportional odds and a non-proportional
# odds population model.  Check how 3 kinds of residuals detect
# non-prop. odds
set.seed(71)
n <- 400
x <- rnorm(n)

par(mfrow=c(2,3))
for(j in 1:2) {     # 1: prop.odds   2: non-prop. odds
  if(j==1) 
    L <- matrix(c(1.4,.4,-.1,-.5,-.9),nrow=n,ncol=5,byrow=TRUE) + x/2 else {
	  # Slopes and intercepts for cutoffs of 1:5 :
	  slopes <- c(.7,.5,.3,.3,0)
	  ints   <- c(2.5,1.2,0,-1.2,-2.5)
      L <- matrix(ints,nrow=n,ncol=5,byrow=TRUE)+
           matrix(slopes,nrow=n,ncol=5,byrow=TRUE)*x
	}
  p <- plogis(L)
  # Cell probabilities
  p <- cbind(1-p[,1],p[,1]-p[,2],p[,2]-p[,3],p[,3]-p[,4],p[,4]-p[,5],p[,5])
  # Cumulative probabilities from left to right
  cp  <- matrix(cumsum(t(p)) - rep(0:(n-1), rep(6,n)), byrow=TRUE, ncol=6)
  y   <- (cp < runif(n)) \%*\% rep(1,6)


  f <- lrm(y ~ x, x=TRUE, y=TRUE)
  for(cutoff in 1:5)print(lrm(y>=cutoff ~ x)$coef)


  print(resid(f,'gof'))
  resid(f, 'score', pl=TRUE)
  # Note that full ordinal model score residuals exhibit a
  # U-shaped pattern even under prop. odds
  ti <- if(j==2) 'Non-Proportional Odds\nSlopes=.7 .5 .3 .3 0' else
    'True Proportional Odds\nOrdinal Model Score Residuals'
  title(ti)
  resid(f, 'score.binary', pl=TRUE)
  if(j==1) ti <- 'True Proportional Odds\nBinary Score Residuals'
  title(ti)
  resid(f, 'partial', pl=TRUE)
  if(j==1) ti <- 'True Proportional Odds\nPartial Residuals'
  title(ti)
}
par(mfrow=c(1,1))


\dontrun{
# Get data used in Hosmer et al. paper and reproduce their calculations
v <- Cs(id, low, age, lwt, race, smoke, ptl, ht, ui, ftv, bwt)
d <- read.table("http://www.umass.edu/statdata/statdata/data/lowbwt.dat",
                skip=6, col.names=v)
d <- upData(d, race=factor(race,1:3,c('white','black','other')))
f <- lrm(low ~ age + lwt + race + smoke, data=d, x=TRUE,y=TRUE)
f
resid(f, 'gof')
# Their Table 7 Line 2 found sum of squared errors=36.91, expected
# value under H0=36.45, variance=.065, P=.071
# We got 36.90, 36.45, SD=.26055 (var=.068), P=.085
# Note that two logistic regression coefficients differed a bit
# from their Table 1
}
}
\keyword{models}
\keyword{regression}
\concept{logistic regression model}
\concept{model validation}
