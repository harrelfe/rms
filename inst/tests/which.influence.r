## Fromm Yuwei Zhu
require(rms)
n <- 100     # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
outcome        <- sample(c(1,2,3), n,TRUE)
sex            <- factor(sample(c('female','male'), n,TRUE))
ddist <- datadist(outcome, age,  blood.pressure, sex)
options( datadist = 'ddist')

modelform <- outcome ~ age+sex+blood.pressure
f <-lrm(modelform, x = TRUE, y = TRUE)
which.influence(f) 
