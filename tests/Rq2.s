## Check ability of quantile regression to estimate stratified medians
require(quantreg)
sm <- function(n) {
  y <- exp(rnorm(n))
  x <- c(rep(0, n/2), rep(1, n/2))
  y[x==1] <- y[x==1] * 1.4
  f <- rq(y ~ x)
  sampmed <- tapply(y, x, median)
  qrmed   <- c(coef(f)[1], sum(coef(f)))
  print(rbind(sampmed,qrmed))
}
for(i in 1:10) sm(100)
for(i in 1:10) sm(1000)

## Compare standard error of mean | x=0 with standard error from quantile
## regression using 4 methods

cse <- function(n) {
  y <- rnorm(n)
  x <- c(rep(0, n/2), rep(1, n/2))
  sem <- sd(y[x==0])/sqrt(n/2)
  semr <- sem*sqrt(pi/2)
  print(c(n=n, sem=sem, semr=semr)) # semr: asympt. se(median)
  f <- rq(y ~ x)
  for(m in c('iid', 'nid', 'ker', 'boot')) { # nid is default
    cat(m, '\n')
    s <- coef(summary(f, se=m))['(Intercept)','Std. Error']
    cat('se method:', m, '  se=', s, '\n')
  }
}
for(i in 1:10) cse(100)
for(i in 1:10) cse(5000)
# nid does appear to work best
