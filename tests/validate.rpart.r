require(rms)
require(rpart)
n <- 100
set.seed(1)
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
y  <- 1*(x1+x2+rnorm(n) > 1)
table(y)
f <- rpart(y ~ x1 + x2 + x3, model=TRUE)
v <- validate(f)
v    # note the poor validation
par(mfrow=c(1,2))
plot(v, legendloc=c(.2,.5))
par(mfrow=c(1,1))


## From http://stackoverflow.com/questions/37053654   Rohan Adur
set.seed(4)
dat = data.frame(X1 = sample(x = c(1,2,3,4,5), size = 100, replace=TRUE))
dat$t = rexp(100, rate=dat$X1)
dat$t = dat$t / max(dat$t)
dat$e = rbinom(n = 100, size = 1, prob = 1-dat$t )

f = rpart(Surv(t, event = e) ~ X1 , data = dat, model=TRUE,
          control=rpart.control(minsplit=30, cp=0.01))
plot(f); text(f)
v <- validate(f)
v
plot(v, legendloc=c(.6,.2))
