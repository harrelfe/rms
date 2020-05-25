if(! require(mice)) quit(save='no')
require(rms)
set.seed(1)
n <- 50
d <- data.frame(x1=runif(n), x2=sample(c('a','b','c'), n, TRUE),
                x3=sample(c('A','B','C','D'), n, TRUE),
                x4=sample(0:1, n, TRUE),
                y=runif(n))
d$x1[1:5]  <- NA
d$x2[3:9]  <- NA
d$x3[7:14] <- NA

a <- aregImpute(~ x1 + x2 + x3 + x4 + y, data=d)
ols(y ~ x1 + x2 + x3 + x4, data=d)

fit.mult.impute(y ~ x1 + x2 + x3 + x4, ols, a, data=d)  # works

m <- mice(d)
d1 <- complete(m, 1)
## ols(y ~ x1 + x2 + x3 + x4, data=d1)  # fails

w <- d1
attr(w$x2, 'contrasts') <- NULL
attr(w$x3, 'contrasts') <- NULL
ols(y ~ x1 + x2 + x3 + x4, data=w)   # works

fit.mult.impute(y ~ x1 + x2 + x3 + x4, ols, m, data=d)

