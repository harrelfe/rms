require(rms)
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
label(age) <- "Age"
sex <- factor(sample(c('male','female'), n, TRUE))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='female'))
dt <- -log(runif(n))/h
label(dt) <- 'Follow-up Time'
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)
units(dt) <- "Year"
sex2 <- factor(as.vector(sex),levels=c('male','female'))
dd <- datadist(age, sex, sex2)
options(datadist='dd')
S <- Surv(dt,e)

f <- npsurv(S ~ sex)
survplot(f, n.risk=TRUE)

f2 <- npsurv(S ~ sex2)
survplot(f2, n.risk=TRUE)

f <- cph(S ~ strat(sex2), surv=TRUE)
survplot(f, n.risk=TRUE, conf.int=.95)

f <- cph(S ~ sex2, surv=TRUE)
survplot(f, n.risk=TRUE, conf.int=.95)
