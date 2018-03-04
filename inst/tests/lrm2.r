# https://github.com/harrelfe/rms/issues/55

require(rms)
set.seed(1)
n <- 20
X2 <- factor(c(rep(0, n/2), rep(1, n/2)))
X21 <- rep(1 : (n/2), 2)
y <- rep(0 : 1, n/2)
options(rmsdebug=TRUE)
f <- lrm(y ~ X2 + X21, method='model.frame')
attributes(f)$Design$mmcolnames

# Problem is inherent to R
colnames(model.matrix(~ X2 + X21))
