require(rms)
set.seed(1)
n <- 100
y <- 0:99
x <- (y / 100) + runif(n)
i <- order(x)
x <- x[i]
y <- y[i]
f <- orm(y ~ x, x=TRUE, y=TRUE)
f
survest(f, data.frame(x=1), times=50)
g     <- orm.fit(x, y)
iref  <- f$interceptRef
k     <- num.intercepts(f)
c(k, iref)
alpha <- coef(f)[1 : k]
beta  <- coef(f)[-(1 : k)]
c(alpha[iref], beta)
lp    <- alpha[iref] + beta * x
flp   <- f$linear.predictors
range(lp - flp)
plogis(alpha['y>=51'] + beta)
lp    <- alpha[iref] + beta
survest(f, linear.predictors=lp, times=50, conf.int=0)
survest(g, linear.predictors=lp, times=50, conf.int=0)
Sf <- Survival(f)
Sg <- Survival(g)   # identical except for covariate name x vs x[1]
# Stores intercepts and betas from initial model fit
# survest runs Survival() from this fit
Sf(50, lp)
Sg(50, lp)

options(calibrate.debug=FALSE)
cal <- calibrate(f, u=50, B=100)
#                 val.surv.args=list(method='smoothkm'),
#                 pred=c(.2, .5, .8)) #seq(0, 1, by=0.1))


