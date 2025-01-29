require(rms)
set.seed(1)
n <- 200
x    <- runif(n)
cens <- round(runif(n, 0.5, 3), 2)
etime <- round(rexp(n), 2)
range(etime[etime <= cens])
# xless(cbind(cens, etime, first=pmin(cens, etime), rt=ifelse(etime < cens, etime, maxt)))
y1 <- pmin(cens, etime)
y2 <- ifelse(etime < cens, etime, Inf)
# source('~/r/rms/r/Ocens.r')
y <- Ocens(y1, y2)
lev <- attr(y, 'levels')
cbind(cens, etime, y1, y2, y, lev[y[,1]], lev[y[,2]])[36, ]
cat('k=', length(lev) - 1, '\n')
f <- orm.fit(y=y, family='logistic', trace=1)
alpha <- unname(coef(f)[1 : num.intercepts(f)])

s <- km.quick(Surv(y1, etime < cens), interval='>=', type='kaplan-meier')
all.equal(f$yunique, s$time)
cdf <- plogis
# cdf <- function(x) exp(-exp(-x))
all.equal(cdf(alpha), s$surv[-1])
range(plogis(alpha) - s$surv[-1])

require(rms)
set.seed(1)
x1 <- rnorm(100)
y  <- x1 + rnorm(100)
ev <- sample(0:1, 100, TRUE)
y[ev == 0] <- pmax(y[ev == 0], 0)
range(y[ev ==1])
Y <- Ocens(y, ifelse(ev, y, Inf))
cbind(y, ev, Y)
f <- orm(Y ~ x1, family='loglog', y=TRUE, lpe=TRUE)
g <- cph(Surv(y, ev) ~ x1)
ordESS(f)

# Add some interval-censored observations
y1 <- y
y2 <- ifelse(ev, y, Inf)
x1 <- c(x1, rnorm(5))
y1 <- c(y1, c(-2, -1.5, -1, -.5, 0))
y2 <- c(y2, c(1, 1, 1, 1, 1))
Y <- Ocens(y1, y2)
f <- orm(Y ~ x1, family='loglog', y=TRUE, lpe=TRUE)
ordESS(f)


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
system.time(km.quick(Surv(y)))                # 0.037s used for initial alphas
system.time(orm.fit(x, y, compstats=FALSE))   # 0.9s

# require(MLEcens)
# w <- cbind(1:N, c(1:10, 14, 12:N), rep(0, N), rep(1, N))
# g <- computeMLE(w, B=c(1,1))   Extremely slow; never finished

# Compare semiparametric model estimates with Turnbull estimates
set.seed(1)
N <- 300
i <- 100 : 200
y <- 1 : N
y2 <- y
y2[i] <- y2[i] + sample(5 : 20, length(i), TRUE)
S <- Surv(y, y2, type='interval2')
# survfit takes 1.3s for N=1000 with 100 interval censored observations
system.time(f <- survival::survfit(S ~ 1, conf.type='none'))
plot(f$time, f$surv, type='s')
# i <- f$time > 90 & f$time < 110
# with(f, plot(time[i], surv[i], type='l'))  # skip betw 99 - 109.5 with y2[i]+10
# abline(v=100)
sum(y == y2)    # gap 100-200
# cbind(y, y2)
Y <- Ocens(y, y2)
# options(orm.fit.debug=3)
system.time(g <- orm.fit(y=Y, compstats=FALSE))  #  0.2s
ti <- g$yunique
s <- c(1, plogis(coef(g)))
lines(ti, s, type='s', col='red')
lines(ti, s, col='blue')
lowest_ic <- min(y[y != y2])
highest_ic <- max(y[y != y2])
c(lowest_ic, highest_ic)
lines(ti[ti < lowest_ic], s[ti < lowest_ic], type='s', col='green')
lines(ti[ti > highest_ic + 1], s[ti > highest_ic + 1], type='s', col='green')
lines(c(lowest_ic - 1, highest_ic + 2), s[ti %in% c(lowest_ic - 1, highest_ic + 2)], col='purple')
