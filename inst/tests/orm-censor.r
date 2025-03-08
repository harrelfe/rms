require(rms)
require(icenReg)
set.seed(1)
n <- 100
x    <- runif(n)
cens <- round(runif(n, 0.5, 3), 2)
etime <- round(rexp(n), 3)
range(etime[etime <= cens])
# xless(cbind(cens, etime, first=pmin(cens, etime), rt=ifelse(etime < cens, etime, maxt)))
y1 <- pmin(cens, etime)
y2 <- ifelse(etime < cens, etime, Inf)
S  <- Surv(y1, etime < cens)
# source('~/r/rms/r/Ocens.r')
y <- Ocens(y1, y2)
# Reconstruct Surv object from Ocens object and check agreement
# with K-M estimates even though positions of right-censored points
# backs up to previous uncensored point
f <- km.quick(S)
S2 <- Ocens2Surv(y)
f2 <- km.quick(S2)
with(f, plot(time, surv, type='s'))
with(f2, lines(time, surv, type='s', col='red'))
m <- min(length(f2$time), length(f$time))
max(abs(f$time[1:m] - f2$time[1:m]))
max(abs(f$surv[1:m] - f2$surv[1:m]))

#lev <- attr(y, 'levels')
# cbind(cens, etime, y1, y2, y, lev[y[,1]], lev[y[,2]])[36, ]
#cat('k=', length(lev) - 1, '\n')
f <- orm.fit(y=y, family='logistic', trace=1)
# Had already concerged at starting values (KM estimates)
# So MLE = K-M
plotIntercepts(f)
alpha <- unname(coef(f)[1 : num.intercepts(f)])

set.seed(2)
x1 <- rnorm(100)
y  <- round(x1 + rnorm(100), 3)
ev <- sample(0:1, 100, TRUE)
y[ev == 0] <- pmax(y[ev == 0], 0)
range(y[ev ==1])

Y <- Ocens(y, ifelse(ev == 1, y, Inf))
cbind(y, ev, Y)
f <- orm(Y ~ x1, family='loglog', y=TRUE, lpe=TRUE)
g <- cph(Surv(y, ev) ~ x1)
ordESS(f)

# S <- cbind(y, ifelse(ev == 1, y, Inf))
# d <- data.frame(S=I(S), x1)
# f <- ic_sp(S ~ x1, data=d, model='po')

# Add some interval-censored observations
y1 <- y
y2 <- ifelse(ev, y, Inf)
range(y1[y1 == y2])
i <- order(y1)
y1 <- y1[i]
y2 <- y2[i]
x1 <- c(x1, rnorm(5))
y1 <- c(y1, c(-2, -1.5, -1, -.5, 0))
y2 <- c(y2, c(1, 1, 1, 1, 1))
# xless(cbind(y1, y2))
f <- Ocens(y1, y2)
# table(attr(f, 'levels'))

Y <- Ocens2ord(Ocens(y1, y2))
f <- attr(Y, 'npsurv')
plot(f$time, f$surv, type='s')   # may be P(T >= t)
S <- Surv(y1, y2, type='interval2')
g <- survfit(S ~ 1)
lines(g$time, g$surv, type='s', col='blue')  #  P(T > t)
h <- Ocens2ord(Ocens(y1, y2))
a <- attributes(h)
np <- a$npsurv
lines(np$time, np$surv, type='s', col='red')

# d <- data.frame(S=I(Surv(y1, y2, type='interval2')), x1)
# f <- ic_sp(S ~ x1, data=d, model='po')
# Bug: falsely claims a y1 is > a y2

d <- data.frame(y1, y2, unc=is.finite(y2), rt=is.infinite(y2))
d$ic <- (1 : 105) > 100
i <- order(d$y1)
d <- d[i, ]
d$i <- 1 : nrow(d)
with(d, plot(i[unc], y1[unc], xlab='', ylab='t'))
with(d, points(i[rt], y1[rt], col='red'))
with(subset(d, ic), for(j in 1:5) lines(c(i[j], i[j]), c(y1[j], y2[j])))

with(d, plot(i[unc], y1[unc], xlab='', ylab='t'))
Y <- Ocens2ord(Ocens(y1, y2))
a <- attributes(Y)
np <- a$npsurv
which(diff(np$surv) >= 0)
# with(a, cbind(levels, upper, ifelse(levels != upper, 1, 0), np$time, np$surv))
# Row 35: somewhat large gap without an intervening uncensored observation
# See https://chatgpt.com/c/679b8d9e-ef84-800a-b216-0ef6ba908a93
# i <- order(y1)
# cbind(y1, y2)[i,]
# f <- orm(Y ~ x1, family='loglog', y=TRUE, lpe=TRUE, trace=1)
# ordESS(f)


# Create a simple dataset without censoring, then add one interval censored value and
# see if information matrix and model chi-square changed very little
set.seed(2)
x <- runif(100)
y <- 1 : 100
f <- orm(y ~ x)
f
y <- Ocens(c(y, 2), c(y, 99))
x <- c(x, runif(1))
g <- orm(y ~ x, lpe=TRUE, y=TRUE)
vcov(f)
vcov(g)
a <- f$info.matrix
b <- g$info.matrix
range(a$ab - b$ab)
range(a$b  - b$b)
ordESS(g)

# For a large dataset add one interval censored observation to check how much longer
# orm.fit runs by using a general sparse matrix for the intercept hessian matrix

set.seed(2)
N <- 100000
x <- matrix(rnorm(N*20), N, 20)
y <- 1 : N
system.time(orm.fit(x, y, compstats=FALSE))   # 0.9s
x <- rbind(x, rnorm(20))
y <- Ocens(c(y, 100), c(y, 150))
system.time(Y <- Ocens2ord(y, verbose=TRUE))  # 0.13
system.time(orm.fit(x, y, compstats=FALSE)) # 1.0s

# require(MLEcens)
# w <- cbind(1:N, c(1:10, 14, 12:N), rep(0, N), rep(1, N))
# g <- computeMLE(w, B=c(1,1))   Extremely slow; never finished

# Compare semiparametric model estimates with Turnbull estimates
set.seed(1)
N <- 40
i <- 10 : 20
y <- 1 : N
y2 <- y
y2[i] <- y2[i] + sample(3 : 10, length(i), TRUE)

y <- Ocens(y, y2)
s <- Ocens2ord(y, nponly=TRUE)  # get initial Turnbull intervals
Y <- Ocens2ord(y)
plot(s$time, s$surv, type='s', xlab='t', ylab='S(t)')

np <- attr(Y, 'npsurv')   # after consolidation
with(np, data.frame(time, surv))
g <- orm.fit(y=y, trace=2, opt_method='LM')  # did not work with 'NR'
ti <- g$yunique
s <- c(1, plogis(coef(g)))
lines(ti, s, type='s', col='red')
