# Thanks: Chris Andrews  <chrisaa@med.umich.edu>
require(survival)
left <- c(1, 3, 5, NA) 
right <-c(2, 3, NA, 4) 
Surv(left, right, type='interval2') 
survreg(Surv(left, right, type='interval2') ~ 1) 

require(rms)
Surv(left, right, type='interval2') # err
args(Surv)
psm(Surv(left, right, type='interval2') ~ 1)
