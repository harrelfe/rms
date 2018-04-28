# From Tamas Ferenci

require(rms)
set.seed(1)
n <- 100
SimData <- data.frame( x1 = sample(c(-1,0,1), n, TRUE),
            x2 = sample(c(-1,0,1), n, TRUE),
            exposure = rlnorm(100) )
SimData$y <- round(runif(n), 2)
dd <- datadist(SimData)
options(datadist="dd")
f <- lrm(y ~ x1 + x2 + offset(log(exposure)), data=SimData, eps=.0001)
d <- orm(y ~ x1 + x2 + offset(log(exposure)), data=SimData, eps=.0001)
max(abs(coef(f) - coef(d)))

h <- orm(y ~ x1 + x2, family='cloglog', data=SimData)
k <- orm(y ~ x1 + x2 + offset( log( exposure ) ),
         family='cloglog', data=SimData)
coef(h) - coef(k)
