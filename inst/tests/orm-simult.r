# From Matthew Shun-Shin <m@shun-shin.com>  2018-01-14

require(rms)
set.seed(1)
m <- 50
d <- expand.grid(arm=c('a','b','c'), i=1 : m)
d$x <- runif(nrow(d))
d$y <- rnorm(nrow(d))
dd <- datadist(d)
options(datadist="dd")

f <- ols(y ~ x + arm, data=d)
summary(f, verbose=TRUE)
summary(f, conf.type='simult', verbose=TRUE)  # simult ignored
#Works
contrast(f, list(arm=c('c','b')), list(arm='a'))
contrast(f, list(arm=c('c','b')), list(arm="a"),
                  conf.type='simultaneous')

g <- orm(y ~ x + arm, data=d)
summary(g, verbose=TRUE)
summary(g, conf.type='simultaneous', verbose=TRUE)  # simult ignored

contrast(g, list(arm=c('b','c')), list(arm='a'))
contrast(g, list(arm=c('b','c')), list(arm='a'), conf.type='simult')

