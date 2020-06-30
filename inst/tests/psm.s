intr <- FALSE  # set to TRUE if running interactivel so xless will run
require(survival)
n <- 400
set.seed(1)
age <- rnorm(n, 50, 12)
sex <- factor(sample(c('Female','Male'),n,TRUE))
# Population hazard function:
h <- .02*exp(.06*(age-50)+.8*(sex=='Female'))
d.time <- -log(runif(n))/h
cens <- 15*runif(n)
death <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)
pol <- function(x, d) cbind(x, x^2)
g <- survreg(Surv(d.time,death) ~ sex*pol(age,2), 
             dist='lognormal')
rg <- residuals(g, type='matrix')[,'dg']

require(rms)
h <- survreg(Surv(d.time,death) ~ sex*pol(age,2), 
             dist='lognormal', x=TRUE)
# lognormal is bad fit for these data
rbind(coef(g), coef(h))

rm(pol)
f <- psm(Surv(d.time,death) ~ sex*pol(age,2), 
              dist='lognormal', x=TRUE, y=TRUE) #, control=survreg.control())
rbind(coef(h), coef(f))
v <- vcov(f, regcoef.only=FALSE)
diag(vcov(h)) / diag(v)

r <- residuals(f, type='matrix')[,'dg']
if(intr) xless(cbind(rg, r))

if(intr) xless(residuals(f, type='score'))
fr <- robcov(f)
diag(vcov(f)) / diag(vcov(fr))


r <- residuals(f)
g <- npsurv(r ~ sex)
survplot(g)

# Generate data where age is irrelevant but PH assumption for sex
# is satisfied (Weibull fits but lognormal doesn't)
set.seed(1)
sex <- factor(sample(c('Female','Male'), n, TRUE))
# Population hazard function:
h <- .02*exp(0.5 + 1.6*(sex=='Female'))
d.time <- -log(runif(n))/h
cens <- 15*runif(n)
death <- ifelse(d.time <= cens,1,0)
d.time <- pmin(d.time, cens)
table(death)

par(mfrow=c(1,2))
for(dist in c('lognormal', 'weibull')) {
f <- psm(Surv(d.time, death) ~ sex,
         dist=dist, x=TRUE, y=TRUE)
r <- residuals(f, type='censored.normalized')
g <- npsurv(r ~ sex)
survplot(g)
lines(r)
title(dist)
}
