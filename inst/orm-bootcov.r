require(rms)
set.seed(1)
n <- 100
y <- sample(1 : 10, n, TRUE)
x1 <- runif(n)
x2 <- runif(n)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001)
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001)
set.seed(1)
fb <- bootcov(f, B=400, eps=.001)
k <- c('y>=5', 'x1', 'x2')
set.seed(1)
gb <- bootcov(g, B=400, eps=.001)
attributes(gb$var)
fb$var[k, k]
gb$var[k, k]

vcov(fb)[k, k]
vcov(gb)
vcov(gb, intercepts='mid')

anova(fb)
anova(gb)
# Still need to understand how bootcov works differently for orm
