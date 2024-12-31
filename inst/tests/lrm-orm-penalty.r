require(rms)

# Example of penalty from help page using lrm:

n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)

f <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
         x=TRUE, y=TRUE)
fu <- f
p <- pentrace(f, seq(.2,1,by=.05))
plot(p)
p$diag      # may learn something about fractional effective d.f. for each original parameter
g <- orm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
         x=TRUE, y=TRUE)
gu <- g

pg <- pentrace(g, seq(.2,1,by=.05))
plot(pg)
cbind(p$diag, pg$diag)

f <- update(fu, penalty=.02)
g <- update(gu, penalty=.02)

cbind(coef(f) / coef(fu), coef(g) / coef(gu))

range(coef(f) - coef(g))
range(vcov(f) - vcov(g))

for(n in c('a', 'b', 'ab')) {
  cat(n, '\n')
  print(f$info.matrix[[n]] / g$info.matrix[[n]])
}

