require(rms)
n <- 400000
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
x4 <- runif(n)
x5 <- runif(n)
x6 <- runif(n)
x7 <- runif(n)
x8 <- runif(n)
x9 <- runif(n)
x10 <- runif(n)
X <- cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
L <- x1 + x2 + x3 - 1.5
y <- ifelse(runif(n) <= plogis(L), 1, 0)
fm <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10
system.time(f <- glm(fm, family=binomial))
print(summary(f), digits=7)
system.time(g <- lrm(fm))
system.time(lrm.fit(X, y))
print(g, digits=7)
coef(f) - coef(g)
sqrt(diag(vcov(f)))/sqrt(diag(vcov(g)))
system.time(h <- orm(fm))
system.time(i <- orm.fit(X, y))
Rprof('orm.fit.out')
of <- orm.fit(X, y)
Rprof(NULL)
system('R CMD Rprof orm.fit.out')

require(MASS)
n <- 300
y <- factor(sample(0:4, n, TRUE))
x1 <- runif(n)
x2 <- runif(n)
x3 <- runif(n)
system.time(f <- polr(y ~ x1 + x2 + x3))
print(summary(f, digits=7))
system.time(g <- lrm(y ~ x1 + x2 + x3))
print(g, digits=7)
c(-f$zeta, f$coefficients) - coef(g)
print( (diag(vcov(f))[c(4:7, 1:3)])/diag(vcov(g)), digits=10)

w <- function(m) {
  x <- runif(200)
  if(m > 0) x[1:m] <- NA
  x
}
set.seed(1)
y <- sample(0:1, 200, TRUE)
x1 <- w(50)
x2 <- w(1)
x3 <- w(2)
x4 <- w(0)
x5 <- w(10)
x6 <- w(11)
x7 <- w(13)
x8 <- w(8)
x9 <- w(7)
x10 <- w(6)
x11 <- w(5)
x12 <- w(4)
x13 <- w(3)
x14 <- w(7)
x15 <- w(18)
x16 <- w(19)
x17 <- w(21)
x18 <- w(23)
x19 <- w(25)
x20 <- w(27)
f <- lrm(y ~ x1+x2+x3+x4+x5+x6+x7+x8+x9+x10+x11+x12+x13+x14+x15+x16+x17+x18+x19+x20)


sink('/tmp/t.tex')
cat('\\documentclass{report}\\usepackage{color,epic,longtable}\\begin{document}',
    sep='\n')
print(f, latex=TRUE)
cat('\\end{document}\n')
sink()

# From Ferenci Tamas  <ferenci.tamas@nik.uni-obuda.hu>
set.seed(1)
d <- data.frame( y = runif( 1000 ) > 0.5, x = rnorm( 1000 ),
                 w = sample( 1:100, 1000, replace = TRUE ) )
wt <- d$w
g <- d[ rep(1 : nrow(d), wt), ]
wtd.mean(d$y, wt); mean(g$y)
wtd.var(d$y, wt); var(g$y)
wtd.mean(d$x, wt); mean(g$x)
wtd.var(d$x, wt); var(g$x)
# The 2 models will disagree if allow to use different knots
k <- quantile(d$x, c(5,20,50,80,95) / 100)
a <- lrm( y ~ rcs( x, k ), data = d, weights = w)
b <- lrm( y ~ rcs( x, k ), data = g )
# xless(a); xless(b)



