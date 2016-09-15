## From Shane McIntosh 

library(rms)

# Prep data
Load(qt50)
dd <- datadist(qt50); options(datadist = "dd")
with(qt50, table(round(entropy, 5)))

# Find smallest model that fails
f <- ols(log1p(post_bugs) ~ entropy, data=qt50, x=T, y=T)
X <- f$x
n <- nrow(X)
y <- f$y
set.seed(1)
# Wild estimates on 22nd resample
for(i in 1 : 22) {
  j <- sample(1 : n, replace=TRUE)
  g <- lm.fit.qr.bare(X[j,], y[j])
  print(coef(g))
}
plot(X[j,], y[j], xlim=range(qt50$entropy))
abline(coef(g))
with(qt50, scat1d(entropy))
# Only 1 value of entropy > .0002, this bootstrap sample had none

set.seed(1)
validate(f, B=100)
