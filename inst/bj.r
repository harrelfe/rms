## From Nicholas Stroustrup nicholas_stroustrup@hms.harvard.edu
require(rms)
deaths = Surv(c(1,2,3,4,5,6),c(1,0,1,1,0,1))
cg = as.factor(as.character(c(1,1,1,0,0,0)))
bj(deaths ~ cg)
