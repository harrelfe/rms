require(rms)
n <- 40
set.seed(1)
y  <- runif(n)
x1 <- runif(n)
x2 <- rnorm(n)
g  <- sample(letters[1:4], n, TRUE)
dd <- datadist(x1, x2, g); options(datadist='dd')

f <- ols(y ~ pol(x1, 2) * g + x2)
pol2 <- function(x) {
  z <- cbind(x, xsq=x^2)
  attr(z, 'nonlinear') <- 2
  z
}
h <- ols(y ~ gTrans(x1, pol2) * g + x2)
specs(h, long=TRUE)
rbind(coef(f), coef(h))
summary(f)
summary(h)
ggplot(Predict(f))
ggplot(Predict(h))

k1 <- list(x1=c(.2, .4), g='b')
k2 <- list(x1=c(.2, .4), g='d')
contrast(f, k1)
contrast(h, k1)
contrast(f, k1, k2)
contrast(h, k1, k2)

f <- ols(y ~ lsp(x1, c(.2, .4)))
# Duplicate lsp but give custom names for columns
lspline <- function(x) {
  z <- cbind(x, x.2=pmax(x - .2, 0), x.4=pmax(x - .4, 0))
  attr(z, 'nonlinear') <- 2:3
  z
}
h <- ols(y ~ gTrans(x1, lspline))
rbind(coef(f), coef(h))
ggplot(Predict(f))
ggplot(Predict(h))
anova(f)
anova(h)

yl <- c(-0.25, 1.25)

# Fit a straight line from x1=0.1 on, but force a flat relationship for x1 in [0, 0.1]
# First do it forcing continuity at x1=0.1
h <- ols(y ~ pmax(x1, 0.1))
xseq <- c(0, 0.099, 1, 0.101, seq(0.2, .8, by=0.1))
ggplot(Predict(h, x1=xseq))
# Now allow discontinuity without a slope change
flin <- function(x) cbind(x < 0.1, x)
h <- ols(y ~ gTrans(x1, flin))
ggplot(Predict(h, x1=xseq), ylim=yl) + geom_point(aes(x=x1, y=y), data=data.frame(x1, y))
# Now have a discontinuity with a slope change
flin <- function(x) cbind(x < 0.1, pmax(x - 0.1, 0))
h <- ols(y ~ gTrans(x1, flin))
ggplot(Predict(h, x1=xseq), ylim=yl) + geom_point(aes(x=x1, y=y), data=data.frame(x1, y))

# Discontinuous linear spline
dlsp <- function(x) {
  z <- cbind(x, x >= 0.2, pmax(x - .2, 0), pmax(x - .4, 0))
  attr(z, 'nonlinear') <- 2:4
  z
}
h <- ols(y ~ gTrans(x1, dlsp))
ggplot(Predict(h), ylim=yl)
ggplot(Predict(h, x1=c(.1, .199, .2, .201, .3, .4, 1)), ylim=yl)

dlsp <- function(x) {
  z <- cbind(x, x >= 0.2, pmax(pmin(x, 0.6) - .2, 0), pmax(pmin(x, 0.6) - .4, 0))
  attr(z, 'nonlinear') <- 2:4
  z
}
h <- ols(y ~ gTrans(x1, dlsp))
ggplot(Predict(h), ylim=yl)

# Try on a categorical predictor
gr <- function(x) cbind(bc=x %in% c('b','c'), d=x == 'd')
h <- ols(y ~ gTrans(g, gr))
ggplot(Predict(h, g))
