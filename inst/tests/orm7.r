require(rms)

x <- 1 : 10
y <- c(0, 1, 0, 0, 0, 1, 0, 1, 1, 1)
lrm(y ~ x)
orm(y ~ x)
orm(y ~ x, family='probit')
orm(y ~ x, family='loglog', trace=1)

x <- 1 : 10
y <- c(0, 2, 0, 1, 0, 2, 2, 1, 1, 2)
f <- orm(y ~ x)
g <- lrm(y ~ x)
coef(f) - g$coef
f$info.matrix
g$info

x <- 1:10
y <- c(0, 2, 0, 1, 0, 3, 2, 1, 1, 3)
f <- orm(y ~ x)
g <- lrm(y ~ x)
orm(y ~ x, family='probit')
require(MASS)
yf <- factor(y)
summary(polr(yf ~ x, method='probit'))
f <- orm(y ~ x, family='loglog')
f$deviance
h <- polr(yf ~ x, method='cloglog')
h
summary(h)


f$info.matrix
g$info.matrix
f$var - g$var[c(1,4), c(1,4)]
coef(f) - coef(g)

n <- 1000
p <- 10
set.seed(1)
x <- matrix(rnorm(n * p), n, p)
y <- 0:999
k <- 999
f <- orm(y ~ x)
g <- lrm(y ~ x)
range(coef(f) - coef(g))
range(vcov(f, intercepts='all') - vcov(g))

require(rms)
n <- 2000; p <- 5; k <- 10
set.seed(1)
x <- matrix(rnorm(n*p), n, p)
y <- sample(0:k, n, TRUE)
# options(rmsdebug=TRUE)
# Get fit object from rms 6.9-0 and makes sure current methods can
# operate on it
f <- readRDS('~/tmp/orm-old-fit.rds')
g <- orm(y ~ x)
class(f$info.matrix); class(g$info.matrix)
f; g
vf <- vcov(f, intercepts='mid'); dim(vf)
vg <- vcov(g, intercepts='mid'); dim(vg)
range(vf - vg)

h <- orm(y ~ x, scale=TRUE)
range(coef(g) - coef(h))
range(vcov(g) - vcov(h))
range(vcov(g, intercepts='all') - vcov(h, intercepts='all'))

g <- update(g, x=TRUE, y=TRUE)
b <- bootcov(g, B=150)
range(vcov(b, intercepts='all') - vcov(g, intercepts='all'))

range(vcov(g, intercepts='all') - b$var)
diag(b$var) / diag(vcov(g, intercepts='all'))

v <- validate(g, B=100)
