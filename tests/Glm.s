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
