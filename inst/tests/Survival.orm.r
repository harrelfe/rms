require(rms)
require(survival)
x <- c(1, 3, 2, 4, 0)
dd <- datadist(x); options(datadist='dd')
y  <- 1:5
y2 <- c(1:4, Inf)
units(y) <- 'day'
Y <- Ocens(y, y2)
f <- orm.fit(y=Y, trace=1, family='loglog')
s <- Survival(f)
s(conf.int=.95)
ti <- c(1, 2, 2.1, 2.2, 2.99999, 3, 4, 4.1)
s(times=ti)
s(times=ti, forcedf=TRUE)
f <- orm(Y ~ x, family='loglog')
s <- Survival(f)
lp <- predict(f, data.frame(x=1:2))
lp
s(lp=lp)
s(lp=lp, time=3)
survest(f, times=2, conf.int=0)
predict(f, data.frame(x=1:2), type='fitted')
# For Y>=2
coef(f)[1] + coef(f)['x'] * (1:2)
s(times=c(1, 3), lp=lp)

X <- predict(f, data.frame(x=1:2), type='x')
s(X=X)
s(conf.int=0.95, X=X)
survest(f, x=X)
survplot(f, x=1:2, conf.int=.95)

f <- survfit(Surv(y, y < 5) ~ 1, conf.type='log-log')
summary(f)

y <- 1:10
x <- rep(0:1, 5)
f <- orm(y ~ x)
f
kint <- f$interceptRef
kint
S <- Survival(f)
alpha <- coef(f)[1 : num.intercepts(f)]
beta  <- coef(f)['x']
alpha
alpha_ref <- alpha[kint]
alpha_ref
lp <- alpha_ref + beta * x
rbind(lp, f$linear.predictors)
lp <- alpha_ref + beta * 1.0
# S(6 | x = 1)
plogis(lp - alpha_ref + alpha['y>=7'])
S(6, lp)
survest(f, data.frame(x=1), times=6)
# S(6 | all original x)
plogis(alpha_ref + beta * x - alpha_ref + alpha['y>=7'])
S(6, f$linear.predictors)
