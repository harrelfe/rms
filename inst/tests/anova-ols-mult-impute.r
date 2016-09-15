require(rms)
set.seed(1)
x1 <- rnorm(100)
x2 <- rnorm(100)
y <- x1 + x2 + rnorm(100)
x1[1:10] <- NA
a <- aregImpute(~ y + x1 + x2)
f <- fit.mult.impute(y ~ x1 + x2, ols, a, data=data.frame(x1,x2,y),
                     n.impute=3, fit.reps=TRUE)
## Show how fit.mult.impute estimates sigma^2
s <- 0
for(i in 1 : 3) s <- s + f$fits[[i]]$stats['Sigma']
c(s / 3, f$stats['Sigma'])

anova(f, test='Chisq')
## Make sure the chi-squares and sums of squares are not from one of the models
for(i in 1 : 3) print(anova(f$fits[[i]], test='Chisq'))


