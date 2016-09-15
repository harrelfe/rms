## Compare log-log ordinal model fit with continuation ratio model
## See Biometrika 72:206-7, 1985
set.seed(171)
type <- 1
n <- 300
y <- sample(0:50, n, rep=TRUE)
sex <- factor(sample(c("f","m"), n, rep=TRUE))
age <- runif(n, 20, 50)
sexo <- sex; ageo <- age
require(rms)
f <- orm(y ~ age + sex, family=loglog)
g <- orm(y ~ age + sex, family=cloglog)
h <- orm(-y ~ age + sex, family=loglog)
i <- orm(-y ~ age + sex, family=cloglog)
p <- function(fit) coef(fit)[c('age','sex=m')]
p(f); p(g); p(h); p(i)
for(type in 1:2) {
  u <- cr.setup(if(type==1) y else -y)
  Y <- u$y
  cohort <- u$cohort
  s <- u$subs
  sex <- sexo[s]
  age <- ageo[s]
  j <- lrm(Y ~ cohort + age + sex)
  print(p(j))
}




