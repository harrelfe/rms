# Example agreement of orm's estimated survival curves and icenReg's
# with a series of random datasets with left, right, and interval censoring

require(rms)
require(icenReg)
Y7 <- NULL
midpts <- TRUE
for(irep in 1 : 25) {
cat('\n----------------------------------------------------\nRepetition', irep, '\n')
set.seed(20 + irep)
N <- 250  # was 80
y  <- sample(10 : 50, N, TRUE)
maxy <- max(y)
y2 <- y
# Left censoring
y[1:30] <- -Inf
# Right censoring
y2[31:60] <- Inf
# Interval censoring
y2[61:90] <- pmin(maxy, y[61:90] + sample(2 : 10, 30, TRUE))
# if(irep == 1) saveRDS(list(y, y2), 'orm-censor3-problem-data.rds')
if(irep == 7) {
  i <- order(y)
  Y7 <- cbind(y, y2)[i, ]
  print(Y7)
}

d <- data.frame(y, y2, grp=factor(rep(1, N)))
f <- ic_np(cbind(y, y2) ~ grp, data=d)
cu <- f$scurves[[1]]
su <- cu$S_curves$baseline
tbull <- cu$Tbull_ints
ti <- (tbull[, 1] + tbull[, 2]) / 2
plot(ti, su, type='s', xlab='Time', ylab='Survival', main=paste('Example', irep))
# print(data.frame(s$Tbull_ints, s$S_curves$baseline)[1:5,])
# plot(f, plot_legend=FALSE, main=paste('Example', irep))

s <- Ocens2ord(Ocens(y, y2), nponly=TRUE)
# print(data.frame(time=c(s$time[1], s$time), surv=c(s$surv, NA))[1:5, ])
lines(c(s$time[1], s$time), c(s$surv, NA), type='s', col='green')
Y <- Ocens2ord(Ocens(y, y2), verbose=FALSE)
s <- attr(Y, 'npsurv')
# print(data.frame(time=s$time, surv=s$surv)[1:5,])

lines(c(s$time[1], s$time), c(s$surv, NA), type='s', col='yellow')

options(orm.fit.debug=FALSE)
f <- orm.fit(y=Ocens(y, y2), trace=1)
ti <- if(midpts || ! length(f$yupper)) f$yunique else (f$yunique + f$yupper) / 2
ti <- c(ti[1], ti)
su  <- c(1, plogis(coef(f)), NA)
lines(ti, su, type='s', col='blue', lwd=2)
}

# See what happens when interval consolidation is not done
if(FALSE) {
y  <- Y7[, 1]
y2 <- Y7[, 2]
Y <- Ocens(y, y2)
Y <- Ocens(y, y2, cons='none')
np <- attr(Y, 'npsurv')
with(np, cbind(time, surv, c(NA, diff(surv) == 0)))
# Y7[y %in% 26:27,]
# Compute initial values forcing them to be in order
init <- qlogis(np$surv[-1])
init[12] <- 0.72
f <- orm.fit(y=Y, initial=init, trace=2, maxit=1)
im <- f$info.matrix
a <- im$a
with(a, plot(row, col))
#Y[11:15,]
require(Matrix)
v <- infoMxop(im)
diag(v)
dg <- a$row == a$col
plot(a$row[dg], a$a[dg])
plot(a$row[! dg], a$a[! dg])
d <- data.frame(row=a$row, col=a$col, a=a$a)
subset(d, abs(a) < 0.1)
subset(d, row %in% 10:13 | col %in% 10:13)
vi <- solve(v)
diag(vi)
vi <- as.matrix(vi)
plot(1 : nrow(vi), diag(vi))
co <- cov2cor(vi)
plotCorrM(co)
co[29:32,29:32]   # r=0.994 for (30, 31)
}
