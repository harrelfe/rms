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
require(rms)
h <- survreg(Surv(d.time,death) ~ sex*pol(age,2), 
             dist='lognormal')
# lognormal is bad fit for these data
rbind(coef(g), coef(h))

rm(pol)
f <- psm(Surv(d.time,death) ~ sex*pol(age,2), 
              dist='lognormal') #, control=survreg.control())
rbind(coef(h), coef(f))
