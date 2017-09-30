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


## From Ferenci Tamas <ferenci.tamas@nik.uni-obuda.hu>
set.seed( 1 )
d <- data.frame(x1 = rnorm( 1000 ),
                x2 = sample( 1:2, 1000, replace = TRUE ),
                x3 = sample( 1:2, 1000, replace = TRUE ) )
d <- transform( d, y = x3*( x1 + x2 )+rnorm( 1000 ) )

dd <- datadist( d )
options( datadist = "dd" )

fit <- ols( y ~ x3*( x1 + x2 ), data = d )

p1 <- Predict( fit, x1, x3 )
p2 <- Predict( fit, x2, x3 )
p  <- rbind(x1=p1, x2=p2)
plot(p, groups='x3', varypred=TRUE)

#Now, if you run plot( p1 ) or plot( p2 ), everything is fine. However, in 
#the last call above the panel for the continuous predictor, x1 is
#fine, the same as plot( p1 ), but for the categorical predictor, it is
#something completely different (and wrong, quite fundamentally: the
#two groups do not even appear). 


