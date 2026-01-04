# See https://github.com/harrelfe/rms/issues/167

library(rms)
set.seed(1)
n <- 100
d <- data.frame(y=factor(sample(c('A','B','C'), size=n, replace=T)), x=rnorm(n))

# lrm fit

f0 <- lrm(y ~ x, data=d, x=T, y=T)
f <- robcov(f0)

Predict(f, x=0:1)
Predict(f, x=0:1, kint=1)
Predict(f, x=0:1, kint=2)

# orm fit

g0 <- orm(y ~ x, data=d, x=T, y=T)
g <- robcov(g0)

Predict(g, x=0:1)
Predict(g, x=0:1, kint=1)
Predict(g, x=0:1, kint=2)

