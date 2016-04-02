## Check ability of quantile regression to estimate stratified medians
require(quantreg)
sm <- function(n, eps=1e-6) {
  y <- exp(rnorm(n))
  x <- c(rep(0, n/2), rep(1, n/2))
  y[x==1] <- y[x==1] * 1.4
  qrmed <- matrix(NA, nrow=2, ncol=2,
                  dimnames=list(c('br','fn'),c('x=0','x=1')))
  for(m in c('br','fn')) {
    f <- if(m == 'br') rq(y ~ x, method=m) else rq(y ~ x, method=m, eps=eps)
    qrmed[m,] <- c(coef(f)[1], sum(coef(f)))
  }
  sampmed <- tapply(y, x, median)
  print(rbind(qrmed,sample=sampmed))
}
for(i in 1:10) sm(100)
for(i in 1:10) sm(1000)

## Compare standard err. of mean | x=0 with standard err. from quantile
## regression using 4 methods

cse <- function(n) {
  y <- rnorm(n)
  x <- c(rep(0, n/2), rep(1, n/2))
  sem <- sd(y[x==0])/sqrt(n/2)
  semr <- sem*sqrt(pi/2)
  res <- vector('numeric', 6)
  names(res) <- c('SEMean','Asympt SEMedian','iid','nid','ker','boot')
  res[1:2] <- c(sem, semr)
  f <- rq(y ~ x)
  for(m in c('iid', 'nid', 'ker', 'boot')) { # nid is default
    s <- coef(summary(f, se=m))['(Intercept)','Std. Error']
    res[m] <- s
  }
  print(t(t(round(res,3))))
}
for(i in 1:10) cse(100)
for(i in 1:10) cse(5000)
# nid does appear to work best

## Compare mean squared err. of quantile estimator of median y | x=E
## in 5-sample problem with orm logistic family estimator.  Also include sample quantile
cmse <- function(n) {   # n = # obs per each of 5 samples
  x <- factor(rep(c('a','b','c','d','e'), n))
  y <- rnorm(5*n)
  s <- x == 'e'
  y[s] <- y[s] + 3
  sampmed <- median(y[s])
  f <- rq(y ~ x)
  qrmed <- coef(f)[1] + coef(f)['xe']
  f <- orm(y ~ x, family=probit)
  if(f$fail) return(c(NA, NA, NA))
  qu <- Quantile(f)
  iref <- f$interceptRef
  ormmed <- qu(.5, z <- coef(f)[iref] + coef(f)['x=e'])
  ormmean <- Mean(f)(z)
  c(sampmed=sampmed, qrmed=qrmed, ormmed=ormmed, ormmean=ormmean)
}
require(rms)
mse <- c(0, 0, 0, 0)
n <- 50
B <- 1000
m <- 0
for(i in 1:B) {
  cat(i, '\r')
  ms <- cmse(n)
  if(!is.na(ms[1])) {
    m <- m + 1
    mse <- mse + (ms - 3) ^ 2
  }
}
m
sqrt(mse/m)   # .123 .124 .126 logistic n=100
              # .173 .176 .172 probit n=50
              # .169 .171 .165 .139 probit n=50   .139=rmse for mean from orm
