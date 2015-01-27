require(rms)
set.seed(1)
n <- 100
y <- sample(1 : 10, n, TRUE)
x1 <- runif(n)
x2 <- runif(n)
f <- lrm(y ~ x1 + x2, x=TRUE, y=TRUE)
g <- orm(y ~ x1 + x2, x=TRUE, y=TRUE)
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
# pseudo.dep, partial .0003
# gof .017


