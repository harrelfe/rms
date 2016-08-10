require(rms)
S <- Surv(1:100)
g <- rep('a', 100)
g[seq(1, 100, by=5)] <- 'b'
f <- npsurv(S ~ g)

survplot(f)
survplot(f, conf='diffbands')
survdiffplot(f)


source('~/R/rms/R/survplotp.npsurv.s')
survplotp(f)

