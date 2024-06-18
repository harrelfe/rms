require(survival)
set.seed(1)
d <- data.frame(x=runif(20), dt=runif(20), ev=sample(0:1, 20, TRUE),
                s=sample(c('f', 'm'), 20, TRUE))
f <- coxph(Surv(dt, ev) ~ x + s, data=d)
f$means
s1 <- survfit(f)$surv

require(rms)
g <- cph(Surv(dt, ev) ~ x + s, data=d, x=TRUE, y=TRUE, surv=TRUE)
g$means
with(d, c(mean(x), mean(s == 'm')))   # = g$means, [2] disagrees with f$means
# g$surv.summary agrees with g$surv and survfit.cph(g)$surv
s2 <- g$surv[-1]
s3 <- survest(g)$surv   # s2=s3, neither=s1

u <- data.frame(x=0.7, s='m')
with(survest(g, u), cbind(time, surv))
tt <- d$dt[2]
survest(g, u, time=tt)

with(survfit(f, u), cbind(time, surv))  # = survest(g, u)

