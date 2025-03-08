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
ms <- modelData(formula=y ~ X)
names(ms)
# f <- lrm(y ~ X, x=TRUE, y=TRUE)

X <- pol(x, 2)
ms <- modelData(formula=y ~ X)
names(ms)
sapply(ms, class)

k <- 4
x <- 1:20
d <- data.frame(x)
m <- model.frame(d, formula=~poly(x, 2))
m
names(m)
m[[1]]
modelData(d, formula= ~ rcs(x,k))
d <- list(x=x, k=6)
m <- modelData(d, ~ rcs(x, k))
m
m[[1]]


b <- 0:9
a <- c(0,1, 1, 2, 2, 3, 4, 7, 7, 9)
e <- rep(c(FALSE, TRUE), length=length(a))
Ocens(a, b)
d <- data.frame(a, b, e)
x <- runif(8)
m <- modelData(d, Ocens(a, b) ~ x, subset=1:8)
attributes(m[[1]])
y <- m[[1]]
class(y)
y1 <- y[1:3,]
class(y1)
y1
y1 <- y[1:3, 1]
class(y1)
y1
y1 <- y[,1]
class(y1)
dim(y1)
m <- modelData(d, Surv(a, e) ~ x, subset=1:8)
y <- m[[1]]
class(y)
y[,1]
class(y[,1])

x <- c(rep('a', 10), rep('b', 11), rep('c', 12))
x <- factor(x, c('a', 'b', 'c', 'd'))
table(x)
y <- runif(length(x))
d <- data.frame(x, y)
m <- modelData(d, y ~ x)
attributes(m$x)


## LCAextend package example like this failed
g <- function() {
  d <- data.frame(x=runif(20), y=sample(0:1, 20,TRUE))
  w <- (1:20)/20
  # d$w <- (1:20)/100 will take precedence
  # return(model.frame(y ~ x, weights=as.vector(w), data=d)) # works
  lrm(y ~ x, weights=as.vector(w), data=d)
}
g()
