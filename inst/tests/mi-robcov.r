require(rms)
set.seed(1)
n <- 100
g <- function() ifelse(runif(n) < 0.2, NA, runif(n))
d <- data.frame(id=1:n, x1=g(), x2=g(), x3=g(), y=rnorm(n))
a <- aregImpute(~ x1 + x2 + x3 + y, data=d)
f <- y ~ x1 + x2 + x3
vcov(fit.mult.impute(f, ols, a, data=d))
vcov(fit.mult.impute(f, ols, a, data=d, robust=TRUE))
w <- rbind(d, d, d, d)
a <- aregImpute(~ x1 + x2 + x3 + y, data=w)
vcov(fit.mult.impute(f, ols, a, data=w))
vcov(fit.mult.impute(f, ols, a, data=w, cluster=w$id))
