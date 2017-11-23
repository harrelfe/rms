## Test handling of NAs in original data by robcov
require(rms)
set.seed(1)
n <- 100
x1 <- runif(n)
x2 <- runif(n)
y  <- x1 + x2 + 3 * runif(n)
x1[1] <- NA
x1 <- c(x1, x1)
x2 <- c(x2, x2)
y  <- c(y, y)
clus <- c(1 : n, 1 : n)
f <- ols(y ~ x1 + x2, subset = 1 : n)
vcov(f)
g <- ols(y ~ x1 + x2, x=TRUE, y=TRUE)
vcov(g)
h <- robcov(g, clus)
h
sqrt(diag(h$var))
vcov(h)
vcov(h) / vcov(f)

