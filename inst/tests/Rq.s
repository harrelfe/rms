###########################################################################
### Purpose: Compare predictions for quantile regression between rq() from
###         the quantreg package and wrapper Rq() from the rms package
### Author: Ben Saville
### Date: 7/26/12
###########################################################################

library(quantreg)
library(rms)

### Simulate data
set.seed(1)
y = rnorm(1000,50,5)
age = sample(5:15,size=1000,replace=TRUE)
gender = as.factor(sample(c("male","female"),size=1000,replace=TRUE))
mydat = data.frame(y,age,gender)

#################### Using rq()
## Fit model with rq
k <- attr(rcs(age,4), 'parms')
rq.test = rq(y ~ rcs(age, k) + gender + rcs(age, k)*gender,
  tau=0.50, data=mydat)

## Create dataset for predictions
p.age = rep(5:15,2)
p.gender = as.factor(rep(c("male","female"),each=11))
p.data = data.frame(p.age,p.gender)
names(p.data) = c("age","gender")
      ## Predictions using predict()
rq.preds = cbind(p.data, predict(rq.test, newdata=p.data))
      ## Predictions using X %*% Beta
p.gender.num = as.numeric(p.gender)-1
X.p = cbind(1, rcs(p.age, k), p.gender.num, rcs(p.age, k)*p.gender.num )
rq.preds.XB = X.p %*% rq.test$coefficients

      ## These match!
cbind(rq.preds,rq.preds.XB)

################## Using Rq()
## Fit model with Rq
Rq.test = Rq(y~ rcs(age, k) + gender +  rcs(age, k)*gender,
  tau=0.5, data=mydat)
## prediction using Predict()
Rq.preds = Predict(Rq.test, age=5:15, gender=c("male","female"),conf.int=FALSE)
## Note predict(Rq.test, newdata=p.data)  gives the same values as Predict()
## Using X %*% Beta
Rq.preds.XB = X.p %*% Rq.test$coefficients

## These don't match!
cbind(Rq.preds, Rq.preds.XB)
