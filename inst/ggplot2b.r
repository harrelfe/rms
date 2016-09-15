## From John Woodill: https://github.com/harrelfe/rms/issues/19
require(rms)
dd = data.frame(x1 = 2 + (runif(200) * 6), x12 = 100 + (runif(200) * 6))
dd$y1 = rep(c(1.2, 1.4), each = 100) * dd$x1 + (runif(200) / 5)

ddist <- datadist(dd)
options("datadist" = "ddist")

g <- ols(y1 ~ x1 + x12, data = dd, x = TRUE, y = TRUE)
a <- Predict(g)

h <- ols(y1 ~ I(x1^2) + I(x12^2), data = dd, x = TRUE, y = TRUE)
b <- Predict(h)

p <- rbind(a,b)
s <- ggplot(p, group = ".set.", ggexpr=TRUE)
s
ggplot(p, group=".set.")
