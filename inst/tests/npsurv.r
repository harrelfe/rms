require(rms)
d <- data.frame(time=1:500, death=sample(0:1, 500, TRUE))
f <- npsurv(Surv(time, death) ~ 1, data=d, conf.type='log-log')
g <- function(y) 1 - y
survplot(f, fun=g)
survplot(f, fun=g, conf='bars')
