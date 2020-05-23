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
k <- c('y>=6', 'x1', 'x2')
set.seed(1)
gb <- bootcov(g, B=400, eps=.001)
list(rownames(f$var), rownames(g$var), rownames(fb$var), rownames(gb$var))
attributes(gb$var)
fb$var[k, k]
gb$var[k, k]

vcov(fb)[k, k]
vcov(gb)
vcov(gb, intercepts='mid')

anova(fb)
anova(gb)
# Still need to understand how bootcov works differently for orm


r <- resid(f, 'score') - resid(g, 'score')
apply(r, 2, function(x) max(abs(x)))
fr <- robcov(f)
gr <- robcov(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001, family='loglog')
gr <- robcov(g)
gb <- bootcov(g, B=200)


# Compare against Yuqi Tian's function
source('robcov_Yuqi.r')
vh <- func_robcov(g, cluster=1:100)
vg <- vcov(gr, intercepts='all')

cv <- function(v1, v2) {
  se1 <- sqrt(diag(v1))
  se2 <- sqrt(diag(v2))
  prn(round(se1 / se2, 3))
  prn(max(abs(v1 - v2)))
}
cv(vg, vh)


cs <- function(fit) {
  sc1 <- resid(fit, 'score')
  sc2 <- func_score(fit)
  max(abs(sc1 - sc2))
}
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001)
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001, family='probit')
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001, family='loglog')
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=0.001, family='cauchit')
cs(g)



