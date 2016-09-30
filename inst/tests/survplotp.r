require(rms)
y <- 1 : 100
units(y) <- 'year'
S <- Surv(y)
g <- rep('a', 100)
g[seq(1, 100, by=5)] <- 'b'

tapply(1:100, g, function(x) sum(x >= 50))  # a=41 b=10

tapply(1:100, g, function(x) length(x) / sum(x))
fs <- survfit(S ~ g)
fs
i <- fs$time %in% c(46, 50)
fs$n.risk[i]   # a=41 b=11
z <- qnorm(.975)
sur <- fs$surv[i]
seS <- fs$std.err[i] * sur   # se on S(t) scale instead of log S(t)
with(fs, cbind(n.risk[i], surv[i], lower[i], upper[i], std.err[i], seS, sur - z * seS, sur + z * seS))
# Last 2 columns not supposed to agree with summary.survfit since summary use log S(t)
# as basis for CLs

# summary(f, times=1:100)$time 

s <- summary(fs, times=50)
s    # a=41 b=10

# Manually compute lower and upper half CL at t=50
mean(sur) + c(-1, 1) * 0.5 * z * sqrt(seS[1]^2 +seS[2]^2) # .3774 .6224
# Compare to width of CL for smallest stratum above       # ,2898 .7191


f <- npsurv(S ~ g)

survplot(f)
survplot(f, conf='diffbands')
survdiffplot(f)   # modern art


survplotp(f)
survplotp(f, aehaz=TRUE)
survplotp(f, times=c(50,60))
survplotp(f, aehaz=TRUE, times=c(5,60))

h <- function(y) 1 - y
survplotp(f, fun=h, ylab='Cumulative Incidence')
survplotp(f, fun=h, aehaz=TRUE, times=c(5, 60))
