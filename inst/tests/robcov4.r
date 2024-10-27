require(rms)
set.seed(1)
cl <- sample(letters[1:5], 100, TRUE)
cli <- 1 : 100
x <- matrix(rnorm(200), ncol=2)
y <- rnorm(100)

f <- ols(y ~ x, x=TRUE, y=TRUE)
vcov(robcov(f))
vcov(robcov(f, cli))
vcov(robcov(f, cl))

g <- lm(y ~ x)
require(sandwich)
vcovHC(g, type='HC0')
vcovCL(g, type='HC0', cluster = ~ cli, cadjust=FALSE)
vcovCL(g, type='HC0', cluster = ~ cl,  cadjust=FALSE)
