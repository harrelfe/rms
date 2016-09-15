# From Rob Kushler <kushler@oakland.edu>
require(rms)
require(MASS)
data(VA)  # use VA lung cancer data in MASS package for examples

# add labels to the factors
VA <- within(VA, {
  treat <- factor(treat,labels=c("std","test"))
  cell <- factor(cell,labels=c("Squamous","Small","Adeno","Large"))
  prior <- factor(prior,labels=c("No","Yes")) })
str(VA)

(VAddist <- datadist(VA))
options(datadist="VAddist")

# model for illustration
(VA.cph2 <- cph(Surv(stime,status) ~ treat*(rcs(Karn,4)+cell+prior), VA, x=TRUE, y=TRUE))

par(mfrow=c(1,2))
survplot(VA.cph2,treat,xlim=c(0,400))
title("Karn=60  cell=Small  prior=No")
survplot(VA.cph2,treat,Karn=30,prior="Yes",xlim=c(0,400))
title("Karn=30  cell=Small  prior=Yes")

ss1 <- survest(VA.cph2)
ss1

S <- with(VA, Surv(stime, status))
f <- cph(S ~ (rcs(Karn,3) + prior)^2 + treat*cell, VA)
with(VA, table(treat)); with(VA, table(cell))
f <- cph(S ~ treat*strat(cell), VA)

f <- cph(Surv(stime,status) ~ (treat+rcs(Karn,3)+prior)^2 + treat*strat(cell),
								VA, x=TRUE, y=TRUE)
                                    
