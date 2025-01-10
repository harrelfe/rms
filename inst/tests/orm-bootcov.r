require(rms)
set.seed(1)
n <- 100
yo <- sample(1 : 10, n, TRUE)
table(yo)
y <- ordGroupBoot(yo, aprob=0.9995, B=1000)
table(yo, y)
x1 <- runif(n)
x2 <- runif(n)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE)
set.seed(1)
fb <- bootcov(f, B=400)
set.seed(1)
gb <- bootcov(g, B=400)
range(vcov(fb, intercepts='all') / vcov(gb, intercepts='all'))
list(rownames(vcov(f)), rownames(vcov(g)), rownames(fb$var), rownames(gb$var))
vcov(gb, intercepts='all')
vcov(gb, intercepts='mid')

anova(fb)
anova(gb)
# Still need to understand how bootcov works differently for orm


r <- resid(f, 'score') - resid(g, 'score')
apply(r, 2, function(x) max(abs(x)))
fr <- robcov(f)
gr <- robcov(g)
range(vcov(fr) / vcov(gr, intercepts='all'))

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, family='loglog')
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
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE)
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, family='probit')
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, family='loglog')
cs(g)

g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, family='cauchit')
cs(g)
