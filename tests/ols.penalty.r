## See http://stats.stackexchange.com/questions/104889/k-fold-or-hold-out-cross-validation-for-ridge-regression-using-r/105453?noredirect=1#comment203976_105453

require(rms)
#random population of 200 subjects with 1000 variables 
M <- matrix(rep(0, 200 * 100), 200, 1000)
for (i in 1 : 200) {
  set.seed(i)
  M[i,] <- ifelse(runif(1000) < 0.5, -1, 1)
}
rownames(M) <- 1:200

##random yvars 
set.seed(1234)
u <- rnorm(1000)
g <- as.vector(crossprod(t(M), u))
h2 <- 0.5 
set.seed(234)
y <- g + rnorm(200, mean=0, sd=sqrt((1 - h2) / h2 * var(g)))

myd <- data.frame(y=y, M)
training.id <- sample(1 : nrow(myd), round(nrow(myd) / 2, 0), replace = FALSE)
test.id <- setdiff(1 : nrow(myd), training.id)

myd_train <- myd[training.id,]
myd_test  <- myd[test.id,] 

frm <- as.formula(paste("y~", paste(names(myd_train)[2:100], collapse="+")))
f <- ols(frm, data = myd_train, x=TRUE, y=TRUE)
p <- pentrace(f, seq(.05, 5, by=.05), noaddzero=TRUE)
plot(p)
