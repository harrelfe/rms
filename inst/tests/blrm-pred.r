require(rms)
stanSet()
set.seed(1)
n <- 100
x1 <- rnorm(n)
x2 <- rnorm(n)
y  <- round(x1 + x2 + rnorm(n), .1)
table(y)
dd <- datadist(x1, x2); options(datadist='dd')
f <- orm(y ~ x1 + x2)
n <- data.frame(x1=c(-.5, .3), x2=.75)

g <- blrm(y ~ x1 + x2)
# Just keep the first 8 posterior draws to make later output smaller
# g$draws <- g$draws[1:8,]
stat <- 'median'
predict(f, n)
p <- predict(g, n, posterior.summary=stat)
p
predict(f, n, type='fitted')
p <- predict(g, n, type='fitted', posterior.summary=stat)
p
predict(f, n, type='fitted.ind')
predict(g, n, type='fitted.ind', posterior.summary=stat)
predict(f, n, type='mean')
p <- predict(g, n, type='mean', posterior.summary=stat)
p

pf <- predict(f, n)
pg <- predict(g, n, posterior.summary=stat)
ybar <- Mean(f)
ybar(pf)
ybar <- Mean(g)
ybar(pg)

yq <- Quantile(f)
yq(lp=pf)   # defaults to median
yq <- Quantile(g)
yq(lp=pg)

ep <- ExProb(f)
ep(pf, y=.5)
ep <- ExProb(g)
ep(pg, y=.5)

f <- update(f, x=TRUE, y=TRUE)

a <- seq(-2, 2, by=1)
Predict(f, x1=a)
Predict(g, x1=a)

Predict(f, x1=a, fun='mean', conf.int=0)
Predict(g, x1=a, fun='mean')

ggplot(Predict(f, x1, fun='mean', conf.int=0))
ggplot(Predict(g, x1, fun='mean'))
