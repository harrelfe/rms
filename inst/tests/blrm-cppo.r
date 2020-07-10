require(rms)
stanSet()
set.seed(1)
x <- rnorm(10)
y <- x + rnorm(10)
tt <- function() {
source('~/R/rms/R/blrm.r')
source('~/R/rms/R/rms.s')
source('~/R/rms/R/predictrms.s')
source('~/R/rms/R/stanMisc.r')
source('~/R/rms/R/contrast.s')
source('~/R/rms/R/Predict.s')
}
system.time(f <- blrm(y ~ x, seed=1, iter=10000))   # 2.2 6.2 for fitter 'lrm' 2.5 6.3 for new
# 2.2 6.0 for new with clustering included in code but not clusters in model
coef(orm(y ~ x))
f <- blrm(y ~ x)
f <- blrm(y ~ x, ~x, cppo=function(y) y)

# Analysis of severity of nausea data from Peterson & Harrell (1990)
d0 <- data.frame(tx=0, y=c(rep(0, 43), rep(1, 39), rep(2, 13), rep(3, 22), rep(4, 15), rep(5, 29)))
d1 <- data.frame(tx=1, y=c(rep(0, 7), rep(1, 7), rep(2, 3), rep(3, 12), rep(4, 15), rep(5, 14)))
d <- rbind(d0, d1)
d$tx <- factor(d$tx, 0:1, c('No cisplatin', 'cisplatin'))
dd <- datadist(d); options(datadist='dd')
with(d, table(tx, y))

coef(blrm(y ~ tx, data=d, method='opt'))
g <- function(y) y==5
standata <- blrm(y ~ tx, ~ tx, cppo=g, data=d, standata=TRUE)
f <- blrm(y ~ tx, ~ tx, cppo=g, data=d, method='opt')
# Compute the treatment effect log(OR) for y=1, 2, 3, 4, 5
h <- f$cppo     # normalized version of f
k <- coef(f)
k
# Before intercept correction for tau * zbar:
# 1.030 0.020 -.293 -.980 -1.681 .940 -.250
k[6] + h(1:5) * k[7]   # matches paper's MLE
fp <- blrm(y ~ tx, ~ tx, cppo=g, data=d)   # get posterior samples
rbind(coef(f), coef(fp, 'mode'), coef(fp, 'mean'))
k <- coef(fp)          # posterior means
k[6] + h(1:5) * k[7]   # close to paper
dat <- data.frame(tx=levels(d$tx))
predictrms(fp, dat, second=TRUE)
contrast(fp, list(tx='cisplatin'), list(tx='No cisplatin'), y=1:5)
Predict(fp, tx, ycut=3)
