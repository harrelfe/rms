require(rms)
options(rmsdebug=TRUE)
y1 <- c(1, 3, 5, 7, NA)
y2 <- c(1, Inf, 7, 7, NA)
label(y1) <- 'Time to event'
units(y1) <- 'day'
y <- Ocens(y1, y2)
y
attributes(y)
y[1:4,]
attributes(y[1:4,])

d <- as.data.frame(y)
class(d$y)
dim(d$y)
d$y

d2 <- na.delete(d)
dim(d2$y)
d2$y
class(d2$y)

d$x <- 11:15
d <- modelData(d, formula=y ~ x)
dim(d$y)
d$y
class(d$y)

d2 <- upData(d)
dim(d2$y)
d2$y
class(d2$y)
d2 <- modelData(d2, formula=y ~ x)
class(d2$y)
dim(d2$y)

