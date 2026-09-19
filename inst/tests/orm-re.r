require(rms)
simulate_pairs <- function(J, beta_true, sigma_true, alpha_true) {
  cl  <- factor(rep(1:J, each = 2))
  x   <- rep(c(0, 1), J)
  u   <- rnorm(J, 0, sigma_true)[as.integer(cl)]
  eta <- beta_true * x + u

  K <- length(alpha_true) + 1
  cumprob <- sapply(alpha_true, function(a) plogis(a - eta))
  cumprob <- cbind(cumprob, 1)
  probs   <- t(apply(cbind(0, cumprob), 1, diff))
  y <- apply(probs, 1, function(p) sample(1:K, 1, prob = p))
  y <- factor(y, ordered = TRUE)

  data.frame(y, x, cl)
}

beta_true  <- 1.0
alpha_true <- c(-1, 0.5, 1.8)
J          <- 800

set.seed(1)
dat <- simulate_pairs(J, beta_true, sigma_true = 3.0, alpha_true)
system.time(g <- with(dat, orm.fit(x, y)))
g[.q(iter, u)]
options(agqStep.dump=FALSE)
system.time(f <- with(dat, orm.fit(x, y, cluster = cl, trace=0)))
if(length(f$fail) && f$fail) stop('fit failed')
f[.q(coefficients, sigma, nAGQ, info.matrix, iter, u, fail)]
vc <- infoMxop(f$info.matrix, invert=TRUE)   # or however vcov.orm normally does this
sqrt(diag(vc))   # should give SEs for alpha and beta

require(ordinal)
system.time(h <- clmm2(y ~ x, random=cl, data=dat, nAGQ=15, Hess=TRUE))
summary(h)


n <- 150
set.seed(2)
xx <- runif(n)
yy <- runif(n)
clus <- factor(c(1 : (n/2), 1 : (n/2)))
system.time(orm.fit(xx, yy))
system.time(f <- orm.fit(xx, yy, cluster=clus))   # .34s
a <- f$info.matrix$a
plot(a$row, a$col)
yf <- as.factor(yy)
system.time(g <- clmm2(yf ~ xx, random=clus, nAGQ=15, Hess=TRUE ))  # 28s
w <- f[.q(coefficients, sigma, nAGQ, info.matrix, iter, deviance, u, fail)]
names(w$info.matrix)
last <- function(z) z[length(z)]
data.frame(beta=last(w$coefficients), sigma=w$sigma, iter=w$iter, nAGQ=w$nAGQ, u=max(abs(w$u)), deviance=w$deviance)
data.frame(beta=g$beta, sigma=g$stDev, iter=g$Niter, nAGQ=w$nAGQ, deviance=-2 * g$logLik)

n <- 1000
set.seed(2)
xx <- runif(n)
yy <- runif(n)
clus <- c(1 : (n/2), 1 : (n/2))
system.time(f <- orm.fit(xx, yy, cluster=clus))   # 0.77s
w <- f[.q(coefficients, sigma, nAGQ, info.matrix, iter, deviance, u, fail)]
a <- f$info.matrix$a
length(a$row)  # 3986 vs 500,000
data.frame(beta=last(w$coefficients), sigma=w$sigma, iter=w$iter, nAGQ=w$nAGQ, u=max(abs(w$u)), deviance=last(w$deviance))
h <- orm.fit(xx, yy)
h$deviance

infoMxop(f$info.matrix, i='sigma_parameters')

f <- orm(yy ~ xx + cluster(clus))

# Test mre - mixing weights for sigma1 sigma2

## Simulate ordinal data under orm's exceedance-probability model
##
##   Prob(Y >= y | X, v_i) = plogis(alpha_y + X*beta + w(t)*v_i)
##
## with cutpoints alpha_1 > alpha_2 > ... > alpha_K (Y takes values
## 0, 1, ..., K), a per-cluster standardized random effect
## v_i ~ N(0,1), and a CONTINUOUS random-effect weight
##
##   w(t) = sigma1*(1 - mre(t)) + sigma2*mre(t)
##
## matching orm.fit's sigma1/sigma2 random-effects design: sigma1 is
## the scale at the anchor time t1 (mre(t1)=0 by construction below),
## sigma2 is the scale toward which w(t) continuously blends as t
## moves away from the anchor. sigma1 and sigma2 are fixed constants
## here (the values being simulated FROM), not estimated.

set.seed(1)

## ---- design ----------------------------------------------------
nc <- 1000L                 # number of clusters (subjects)
nt <- 10L                   # observations per cluster
n  <- nc * nt

times <- c(0, 1, 2, 4, 7, 11, 16, 22, 29, 37)  # follow-up times within a
# cluster (t1=0 is the anchor), extended from
# the earlier 5-point schedule in the same
# increasing-gap style; same schedule for
# every cluster here -- jitter or randomize
# per cluster if you want unequal spacing

## ---- true parameter values --------------------------------------
beta   <- 0.6                     # fixed effect of x
alpha  <- c(2, 1, 0, -1, -2)      # K=5 cutpoints, DECREASING (6 ordinal
# categories, Y = 0..5) -- required for
# Prob(Y>=y) to be monotone decreasing in y
sigma1 <- 1.0                     # random-effect scale at the anchor (t=t1)
sigma2 <- 0.4                     # random-effect scale as t moves away from
# the anchor (set negative here instead,
# e.g. sigma2 <- -0.4, to simulate the
# "subtract an over-induced correlation"
# regime instead of the "add" regime above)

## ---- continuous mre(t): 0 at the anchor, smoothly saturating -----
## toward 1 as t grows. This is the continuous analog of the
## recommended discrete "0 at t1, 1 thereafter" step function -- any
## smooth, user-specified function of t with mre(t1)=0 works equally
## well here; this one is a simple saturating exponential ramp.
decay_rate <- 0.5   # larger = faster approach to mre=1 (i.e. to sigma2)
mre_fun <- function(t) 1 - exp(-decay_rate * (t - times[1]))

w_fun <- function(t) sigma1 * (1 - mre_fun(t)) + sigma2 * mre_fun(t)

## ---- assemble the data -------------------------------------------
cluster <- rep(1 : nc, each = nt)
time    <- rep(times,  times = nc)
mre     <- mre_fun(time)
x       <- rnorm(n)                    # a single continuous covariate

v   <- rnorm(nc, 0, 1)                 # per-cluster standardized random effect
v_i <- v[cluster]

eta <- beta * x + w_fun(time) * v_i    # linear predictor

## ---- draw Y directly from the exceedance probabilities -----------
## Pexceed[, j] = Prob(Y >= j | eta) for j = 1, ..., K (columns
## decreasing left to right for any fixed eta, since alpha is
## decreasing). Y = the number of thresholds a single uniform draw
## falls under -- the standard way to sample a cumulative/
## proportional-odds outcome directly from its exceedance
## probabilities, without needing a latent-variable detour.
K <- length(alpha)
Pexceed <- sapply(alpha, function(a) plogis(a + eta))   # n x K matrix
u <- runif(n)
y <- rowSums(u <= Pexceed)                              # values 0..K

dat <- data.frame(cluster = cluster, time = time, mre = mre, x = x, y = y)

## ---- sanity checks -------------------------------------------------
cat("n:", nrow(dat), " clusters:", nc, " obs/cluster:", nt, "\n")
cat("category counts:\n"); print(table(dat$y))
cat("\nmre(t) at the observed times:\n")
print(round(unique(data.frame(time = times, mre = mre_fun(times),
                              w = w_fun(times))), 3))

## dat is now ready to pass to orm.fit()/orm(), e.g.:
##   f <- orm(y ~ x + cluster(cluster), mre = dat$mre, data = dat)
## once your own code for constructing mre from cluster/time is wired
## into orm()'s own formula interface.

set.seed(3)
x2 <- runif(length(x))
f <- orm(y ~ x + x2 + cluster(cluster) + mix_re(mre), trace=1, x=TRUE, y=TRUE)
latex(f)
anova(f)
anova(f, test='LR')
f
dd <- datadist(x, x2); options(datadist='dd')
summary(f, x=c(-2, 2), x2=c(.2, .7))
contrast(f, list(x=2, x2=0), list(x=-2, x2=0))
contrast(f, list(x=2, x2=0), list(x=-2, x2=0), conf.type='profile')
contrast(f, list(x=0, x2=.7), list(x=0, x2=.2))
contrast(f, list(x=0, x2=.7), list(x=0, x2=.2), conf.type='profile')


# Compare deviances of this well-fitting model to an oversimplified random
# effects model and to a model ignoring random effects
g <- orm(y ~ x + x2 + cluster(cluster), x=TRUE, y=TRUE)
g
anova(g)
anova(g, test='LR')

h <- orm(y ~ x + x2, x=TRUE, y=TRUE)
h
anova(h)
anova(h, test='LR')
deviance(f)
deviance(g)
deviance(h)

# 37s
# i <- clmm2(factor(y) ~ x, random=factor(cluster), nAGQ=15, Hess=TRUE)
# -2 * logLik(i)
# deviance(g)  # identical
# summary(i)  # z=31.5720 = orm

