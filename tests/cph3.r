## Explore width of confidence intervals, checking that it is zero
## at the median of the single predictor

require(rms)
set.seed(1)
a <- 1 : 100
dd <- datadist(a); options(datadist='dd')
S <- Surv(a + 100 * runif(100))
f <- cph(S ~ a)
ab <- list(v = median(a))
plot(Predict(f), abline=ab)
f <- cph(S ~ pol(a, 2))
plot(Predict(f), abline=ab)
plot(Predict(f, ref.zero=TRUE), abline=ab)

b <- sample(1:100) + 100
dd <- datadist(a, b)
f <- cph(S ~ pol(a, 2) + b)
plot(Predict(f), abline=ab)
plot(Predict(f, ref.zero=TRUE), abline=list(v=c(median(a),median(b))))

     
