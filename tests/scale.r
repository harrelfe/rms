## Test center option for lrm
require(rms)
x <- rnorm(30)
y <- sample(0:4, 30, TRUE)
f <- lrm(y ~ pol(x, 2), scale=FALSE)
g <- lrm(y ~ pol(x, 2), scale=TRUE)
coef(f) - coef(g)
d <- vcov(f) - vcov(g)
d
max(abs(d))
