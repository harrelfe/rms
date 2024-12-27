# https://discourse.datamethods.org/t/statistically-efficient-ways-to-quantify-added-predictive-value-of-new-measurements/2013/61?u=f2harrell

require(rms)

n <- 1000
set.seed(1234)
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
time           <- runif(n, 1, 15)
status         <- round(age/100,0)
d <- data.frame(age,blood.pressure,cholesterol,time,status)

f <- cph(
        Surv(time, status) ~ rcs(age, 3) * rcs(blood.pressure, 3),
        x = TRUE,
        y = TRUE,
        data = d
    )

g <- bootcov(f, B=50)
rexVar(g, data=d)
