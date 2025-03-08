require(rms)
set.seed(2)
x <- matrix(runif(20*3), 20, 3)
y <- x[, 1] + runif(20)
f <- orm(y ~ x, x=TRUE, y=TRUE)
validate(f, B=100)

cens <- runif(20, 1, 2)
sum(y < cens)
label(y) <- 'Time to event A'
units(y) <- 'day'
Y <- Ocens(y, ifelse(y < cens, y, Inf))
Y
f <- orm(Y ~ x, x=TRUE, y=TRUE)
f
# options(rmsdebug=FALSE, orm.fit.debug=FALSE, orm.fit.debug2=TRUE, validate.debug=TRUE)
validate(f)

