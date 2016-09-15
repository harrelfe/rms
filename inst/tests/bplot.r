require(rms)
x1 <- runif(100)
x2 <- runif(100)
y  <- x1 + 2 * x2 + 3 * runif(100)
dd <- datadist(x1, x2); options(datadist='dd')
f  <- ols(y ~ x1 * x2)
f
p <- Predict(f, x1, x2, np=20)
bplot(p, lfun=wireframe, col='red')
bplot(p, lfun=wireframe, col='red', xlab='Age (days)', xlabrot=-10,
      cex.lab=1.4)
