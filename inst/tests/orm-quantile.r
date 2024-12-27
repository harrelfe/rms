require(rms)
set.seed(1)
n <- 500
x <- runif(n)
y <- round(x + runif(n, -.1, .1), 2)
f <- orm(y ~ x)
f
d <- data.frame(x=c(-.25, 0, .5, 1))
qu <- Quantile(f)
qu(0.5, lp=predict(f, d))
X <- predict(f, d, type='x')
qu(0.5, X=X, conf.int=0.95)
qu(0.75, X=X, conf.int=0.95)

ex <- ExProb(f)
ex(y=0.5, lp=predict(f, d))
ex(y=0.5, X=X, conf.int=0.95)

M <- Mean(f)
M(lp=predict(f, d))
M(X=X, conf.int=0.05)
