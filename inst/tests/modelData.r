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
Xr <- X; class(Xr) <- c('rms', class(Xr))
ms <- rms:::modelData(formula=y ~ X)
names(ms)
# f <- lrm(y ~ X, x=TRUE, y=TRUE)

X <- pol(x, 2)
ms <- rms:::modelData(formula=y ~ X)
names(ms)

k <- 4
x <- 1:20
d <- data.frame(x)
rms:::modelData(d, formula= ~ rcs(x,k))
d <- list(x=x, k=6)
rms:::modelData(d, ~ rcs(x, k))
