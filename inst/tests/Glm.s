require(rms)
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)

f <- Glm(counts ~ outcome + treatment, family=poisson(), x=TRUE, y=TRUE)
g <- bootcov(f,B=100)
f
g
diag(vcov(g))/diag(vcov(f))

x <- runif(1000)
y <- ifelse(runif(1000) < 0.5, 1, 0)
f <- Glm(y ~ x, family=binomial(), x=TRUE, y=TRUE)
g <- bootcov(f, B=100)
diag(vcov(f))/diag(vcov(g))

###########################################################
## Test offset()
## From rfunction.com/archives/223 and Max Gordon
# Setup some variables suited for poisson regression
Y  <- c(15,  7, 36,  4, 16, 12, 41, 15)
N  <- c(4949, 3534, 12210, 344, 6178, 4883, 11256, 7125)
x1 <- c(-0.1, 0, 0.2, 0, 1, 1.1, 1.1, 1)
x2 <- c(2.2, 1.5, 4.5, 7.2, 4.5, 3.2, 9.1, 5.2)
# Setup the rms environment
ddist <- datadist(Y, N, x1, x2)
options(datadist="ddist")

#############################
# Poisson regression #
#############################
form <- Y ~ offset(log(N)) + x1 + x2
a <- Glm(form, family=poisson)
b <- glm(form, family=poisson)
cbind(coef(a), coef(b))

nd <- data.frame(x1=1, x2=1.5, N=c(1, 1000))
cbind(predict(a, nd), predict(b, nd))
Predict(a, x1=1, x2=1.5, offset=list(N=1000))


## Try with lm and ols
a <- ols(form)
b <- lm(form)
cbind(coef(a), coef(b))
cbind(predict(a, nd), predict(b, nd))
Predict(a, x1=1, x2=1.5, offset=list(N=1000))
cbind(fitted(a), fitted(b))
cbind(resid(a), resid(b))
