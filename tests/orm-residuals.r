require(rms)
set.seed(1)
n <- 100
y <- sample(1 : 10, n, TRUE)
x1 <- runif(n)
x2 <- runif(n)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001)
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE, eps=.001)
max(abs(coef(f) - coef(g)))
max(abs(vcov(f) - vcov(g, intercepts='all')))

options(digits=4)
dm <- function(x) if(length(dim(x))) dim(x) else length(x)
for(type in c('li.shepherd', 'ordinary', 'score', 'pearson',
              'deviance', 'pseudo.dep', 'partial',
              'dfbeta', 'dfbetas', 'dffit', 'dffits', 'hat', 'gof', 'lp1')) {
  cat(type)
  rf <- resid(f, type=type)
  cat(' lrm dim', dm(rf))
  rg <- resid(g, type=type)
  cat(' orm dim', dm(rg))
  cat(' max |difference|', max(abs(rf - rg)), '\n')
}

options(digits=7)
diag(vcov(f)) / diag(vcov(g, intercepts='all'))
diag(vcov(robcov(f))) / diag(vcov(robcov(g), intercepts='all'))

rf <- robcov(f)
rg <- robcov(g)
max(abs(rf$var - rg$var))

max(abs(vcov(rf, intercepts='all') - vcov(rg, intercepts='all')))
vcov(rf, regcoef.only=TRUE, intercepts='none')
vcov(rg, regcoef.only=TRUE, intercepts='none')
