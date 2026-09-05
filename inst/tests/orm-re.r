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
system.time(f <- orm.fit(xx, yy, cluster=clus))   # 0.35s
w <- f[.q(coefficients, sigma, nAGQ, info.matrix, iter, deviance, u, fail)]
a <- f$info.matrix$a
length(a$row)  # 3986 vs 500,000
data.frame(beta=last(w$coefficients), sigma=w$sigma, iter=w$iter, nAGQ=w$nAGQ, u=max(abs(w$u)), deviance=last(w$deviance))
h <- orm.fit(xx, yy)
h$deviance

infoMxop(f$info.matrix, i='log(sigma)')

f <- orm(yy ~ xx + cluster(clus))

