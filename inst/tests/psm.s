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
xless(cbind(rg, r))

xless(residuals(f, type='score'))
fr <- robcov(f)
diag(vcov(f)) / diag(vcov(fr))


