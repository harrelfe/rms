\name{residuals.cph}
\alias{residuals.cph}
\title{Residuals for a cph Fit}

\description{
Calculates martingale, deviance, score or Schoenfeld residuals 
(scaled or unscaled) or influence statistics for a
Cox proportional hazards model. This is a slightly modified version
of Therneau's \code{residuals.coxph} function. It assumes that \code{x=TRUE} and
\code{y=TRUE} were specified to \code{cph}, except for martingale
residuals, which are stored with the fit by default.
}
\usage{
\method{residuals}{cph}(object,
      type=c("martingale", "deviance", "score", "schoenfeld", 
             "dfbeta", "dfbetas", "scaledsch", "partial"), \dots)
}
\arguments{
\item{object}{a \code{cph} object}
\item{type}{
character string indicating the type of residual desired;
the default is martingale.
Only enough of the string to determine a unique match is required.
Instead of the usual residuals, \code{type="dfbeta"} may be specified
to obtain approximate leave-out-one \eqn{\Delta \beta}s.  Use
\code{type="dfbetas"} to normalize the \eqn{\Delta \beta}s for
the standard errors of the regression coefficient estimates.
Scaled Schoenfeld residuals (\code{type="scaledsch"}, Grambsch and
Therneau, 1993) better 
reflect the log hazard ratio function than ordinary Schoenfeld
residuals, and they are on the regression coefficient scale.  
The weights use Grambsch and Therneau's "average variance" method.
}
\item{\dots}{see \code{\link[survival]{residuals.coxph}}}
}
\value{
The object returned will be a vector for martingale and deviance 
residuals and matrices for score and schoenfeld residuals, dfbeta, or dfbetas.
There will
be one row of residuals for each row in the input data (without \code{collapse}).
One column of score and Schoenfeld
residuals will be returned for each column in the model.matrix.
The scaled Schoenfeld residuals are used in the \code{cox.zph} function.


The score residuals are each individual's contribution to the score
vector.  Two transformations of this are often more useful: \code{dfbeta} is
the approximate change in the coefficient vector if that observation
were dropped, and \code{dfbetas} is the approximate change in the coefficients,
scaled by the standard error for the coefficients.
}
\references{
T. Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.


P. Grambsch, T. Therneau. "Proportional hazards tests and diagnostics
based on weighted residuals", unpublished manuscript, Feb 1993.
}
\seealso{
\code{\link{cph}}, \code{\link[survival]{coxph}}, \code{\link[survival]{residuals.coxph}}, \code{\link{cox.zph}}, \code{\link{naresid}}
}
\examples{
# fit <- cph(Surv(start, stop, event) ~ (age + surgery)* transplant, 
#            data=jasa1)
# mresid <- resid(fit, collapse=jasa1$id)


# Get unadjusted relationships for several variables
# Pick one variable that's not missing too much, for fit


n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
cens   <- 15*runif(n)
h      <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
d.time <- -log(runif(n))/h
death  <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)


f <- cph(Surv(d.time, death) ~ age + blood.pressure + cholesterol, iter.max=0)
res <- resid(f) # This re-inserts rows for NAs, unlike f$resid
yl <- quantile(res, c(10/length(res),1-10/length(res)), na.rm=TRUE)
# Scale all plots from 10th smallest to 10th largest residual
par(mfrow=c(2,2), oma=c(3,0,3,0))
p <- function(x) {
  s <- !is.na(x+res)
  plot(lowess(x[s], res[s], iter=0), xlab=label(x), ylab="Residual",
       ylim=yl, type="l")
}
p(age); p(blood.pressure); p(cholesterol)
mtext("Smoothed Martingale Residuals", outer=TRUE)


# Assess PH by estimating log relative hazard over time
f <- cph(Surv(d.time,death) ~ age + sex + blood.pressure, x=TRUE, y=TRUE)
r <- resid(f, "scaledsch")
tt <- as.numeric(dimnames(r)[[1]])
par(mfrow=c(3,2))
for(i in 1:3) {
  g <- areg.boot(I(r[,i]) ~ tt, B=20)
  plot(g, boot=FALSE)  # shows bootstrap CIs
}                  # Focus on 3 graphs on right
# Easier approach:
plot(cox.zph(f))    # invokes plot.cox.zph
par(mfrow=c(1,1))
}
\keyword{survival}
\concept{model validation}
