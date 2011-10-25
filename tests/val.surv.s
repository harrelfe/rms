require(rms)
set.seed(123)              # so can reproduce results
n <- 2500
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('Male','Female'), n, rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
t <- -log(runif(n))/h
units(t) <- 'Year'
label(t) <- 'Time to Event'
ev <- ifelse(t <= cens, 1, 0)
t <- pmin(t, cens)
S <- Surv(t, ev)
d <- data.frame(t, ev, sex, age)
 
ddist <- datadist(d); options(datadist="ddist")
 
fit <- cph(Surv(t,ev) ~ sex + age, data=d, subset=1:1500,
           surv=TRUE, x=TRUE, y=TRUE, time.inc=5)
 
vd <- d[-(1:1500),]
vs <- val.surv(fit, vd, S=with(vd, Surv(t, ev)), u=5)
par(mfrow=c(2,1))
plot(vs)
 
g <- survest(fit, vd, times=5)
vs <- val.surv(fit, vd, S=with(vd, Surv(t,ev)), u=5, est.surv=g$surv)
plot(vs)
