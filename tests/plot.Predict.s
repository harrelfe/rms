# From question by Mike Babyak, Duke U
require(rms)
n     = 30 
group = factor(sample(c('a','b','c'), n, TRUE))
x1    = runif(n)
dat   = data.frame(group, x1,
                   y = as.numeric(group) + 0.2*x1 + rnorm(n) )

d <- datadist(dat) ; options(datadist="d") 

f <- ols(y ~ x1 + group, data=dat)
p <- Predict(f, group) 
plot(p, ~group, nlines=TRUE, type='p', ylab='fitted Y', xlab='Treatment',
     pch=4, lwd=3)
p <- Predict(f, x1=seq(0,1,by=.1), group)
plot(p, ~ x1, groups='group', col=3:1)


