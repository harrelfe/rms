## From Vikki

require(rms)

n <- 1000 
set.seed(110222) 
data <- matrix(rep(0, 5000), ncol=5) 
data[, 1] <- sample(1:3, n, rep=TRUE, prob=c(.32,  .30,  .38)) 
for (i in 1:1000) { 
if (data[i, 1] == 1) data[i, 2] <- sample(1:3, 1, prob=c(.76,  .18,  .06)) 
if (data[i, 1] == 2) data[i, 2] <- sample(1:3, 1, prob=c(.67,  .24,  .09)) 
if (data[i, 1] == 3) data[i, 2] <- sample(1:3, 1, prob=c(.47,  .37,  .16))} 
for (i in 1:1000) { 
if (data[i, 1] == 1) data[i, 3] <- sample(1:4, 1, prob=c(.70,  .19,  .03,  .08)) 
if (data[i, 1] == 2) data[i, 3] <- sample(1:4, 1, prob=c(.42,  .28,  .12,  .18)) 
if (data[i, 1] == 3) data[i, 3] <- sample(1:4, 1, prob=c(.11,  .29,  .30,  .30))} 
for (i in 1:1000) { 
if (data[i, 3] == 1) data[i, 4] <- 12*rgamma(1000, rate=0.4, shape=1.7)[c(sample(26:975, 1, prob=c(rep(1/950, 950))))] 
if (data[i, 3] == 2) data[i, 4] <- 12*rgamma(1000, rate=0.9, shape=1.7)[c(sample(26:975, 1, prob=c(rep(1/950, 950))))] 
if (data[i, 3] == 3) data[i, 4] <- 12*rgamma(1000, rate=1.2, shape=0.6)[c(sample(26:975, 1, prob=c(rep(1/950, 950))))] 
if (data[i, 3] == 4) data[i, 4] <- 12*rgamma(1000, rate=1.5, shape=0.7)[c(sample(26:975, 1, prob=c(rep(1/950, 950))))]} 
for (i in 1:1000) { 
if (data[i, 3] == 1) data[i, 5] <- sample(c(0, 1), 1, prob=c(.53,  .47)) 
if (data[i, 3] == 2) data[i, 5] <- sample(c(0, 1), 1, prob=c(.17,  .83)) 
if (data[i, 3] == 3) data[i, 5] <- sample(c(0, 1), 1, prob=c(.05,  .95)) 
if (data[i, 3] == 4) data[i, 5] <- sample(c(0, 1), 1, prob=c(.06,  .94))} 

d <- data.frame(tumor=factor(data[,1]), ecog=factor(data[,2]), rx=factor(data[,3]), os=data[,4], censor=data[,5])
S <- with(d, Surv(os, censor))

## Check collinearity of rx with other predictors
lrm(rx ~ tumor*ecog, data=d)
## What is the marginal strength of rx (assuming PH)?
cph(S ~ rx, data=d)
## What is partial effect of rx (assuming PH)?
anova(cph(S ~ tumor + ecog + rx, data=d))
## What is combined partial effect of tumor and ecog adjusting for rx?
anova(cph(S ~ tumor + ecog + strat(rx), data=d), tumor, ecog) ## nothing but noise
## What is their effect not adjusting for rx
cph(S ~ tumor + ecog, data=d)  ## huge

f <- cph(S ~ tumor + ecog, x=TRUE, y=TRUE, surv=TRUE, data=d)
set.seed(1)
validate(f, B=100, dxy=TRUE)
w <- rep(1, 1000)   #  only one stratum, doesn't change model
## model.matrix no longer works with one stratum
if(FALSE) {
f <- cph(S ~ tumor + ecog + strat(w), x=TRUE, y=TRUE, surv=TRUE, data=d)
set.seed(1)
validate(f, B=100, dxy=TRUE, u=60)
## identical to last validate except for -Dxy
}

f <- cph(S ~ tumor + ecog + strat(rx), x=TRUE, y=TRUE, surv=TRUE, time.inc=60, data=d)
set.seed(1)
validate(f, B=100, u=60)  ## no predictive ability
set.seed(1)
validate(f, B=100, dxy=TRUE, u=60)
## Only Dxy indicates some predictive information; large in abs. value
## than model ignoring rx (0.3842 vs. 0.3177)
