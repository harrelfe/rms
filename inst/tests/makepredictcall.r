if(FALSE) makepredictcall.rms <- function(var, call)
  {
    # rms transformation functions using parms information/argument
    funs <- c('rcs', 'pol', 'lsp', 'catg', 'scored', 'strat', 'gTrans')
    for(f in funs) {
        if(as.character(call)[1L] == f ||
            (is.call(call) && identical(eval(call[[1L]]), get(f)))) {
          cat(f, 'hit\n')
        call <- call[1L:2L]
        call["parms"] <- attributes(var)["parms"]
        break
        }
    }
    call
}

require(rms)
require(survival)
x <- 1:10
set.seed(1)
y <- Surv(runif(10))
dd <- datadist(x); options(datadist='dd')
p <- function(form) {
  nd <- data.frame(x=c(1, 5, 10))
  f <- cph(form, eps=1e-8, iter.max=80)
  t1 <- predict(f, nd, type='terms')
  prn(t1)
  g <- coxph(form, control=coxph.control(iter.max=80))
  prn(attr(g$terms, 'predvars'))
  t2 <- predict(g, nd, type='terms')
  prn(t2)
  prn(t1 - t2)
  prn(predict(f, nd) - predict(g, nd))
}
p(y ~ rcs(x, 4))
p(y ~ lsp(x, 5))
p(y ~ pol(x, 2))

g <- function(x) {
    X <- cbind(x, (x - 5)^2)
    attr(X, 'nonlinear') <- 2
    X
}
p(y ~ gTrans(x, g))
