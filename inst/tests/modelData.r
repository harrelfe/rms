require(rms)


x <- runif(20)
x2 <- runif(20)
X <- cbind(x, x2)
y <- sample(0:1, 20, TRUE)

m <- model.frame(y ~ X)
names(m)
mft <- attr(m, 'terms')
mft
attr(mft, 'term.labels')
p  <- terms(y ~ X)
attr(p, 'term.labels')
all.vars(y ~ X)
d <- data.frame(y, I(X))
names(d)
source('~/R/rms/R/rms.s')
Xr <- X; class(Xr) <- c('rms', class(Xr))
ms <- modelData(formula=y ~ X)
names(ms)
# Makes 3 variables; need 2
f <- lrm(y ~ X, x=TRUE, y=TRUE)

X <- pol(x, 2)
ms <- rms:::modelData(formula=y ~ X)
names(ms)

