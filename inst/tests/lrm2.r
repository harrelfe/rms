# https://github.com/harrelfe/rms/issues/55

require(rms)
set.seed(1)
n <- 100
X2 <- factor(c(rep(0, n/2), rep(1, n/2)))
X21 <- runif(n)
y <- sample(0:1, n, TRUE)
f <- lrm(y ~ X2 + X21)
