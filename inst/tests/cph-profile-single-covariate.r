require(rms)
library(survival)

# Regression test for GitHub issue #173:
# contrast(..., conf.type = 'profile') on single-covariate cph() models
# passes a zero-column design matrix to coxphFit(), where survival's null-model
# fitters omit res$iter.

dd <- datadist(ovarian)
options(datadist = 'dd')

# Helper: exact profile likelihood CI and LR stat for a single-covariate Cox model
# by evaluating survival::coxph() directly at fixed offset theta * x.
exact_profile <- function(formula_full, formula_null_rhs, data, xvec,
                          ties = 'efron', wvec = NULL, conf.int = 0.95) {
  crit <- qchisq(conf.int, 1)
  d0 <- data
  d0$.wt <- if(is.null(wvec)) rep(1, nrow(data)) else wvec
  fit_full <- coxph(formula_full, data = d0, ties = ties, weights = .wt)
  ll_full  <- fit_full$loglik[length(fit_full$loglik)]
  ll_at <- function(theta) {
    d <- d0
    d$.off <- theta * xvec
    f0 <- coxph(update(formula_full, paste('~', formula_null_rhs, '+ offset(.off)')),
                data = d, ties = ties, weights = .wt)
    f0$loglik[length(f0$loglik)]
  }
  lr   <- 2 * (ll_full - ll_at(0))
  est  <- unname(coef(fit_full)[1])
  se   <- sqrt(vcov(fit_full)[1, 1])
  obj  <- function(theta) 2 * (ll_full - ll_at(theta)) - crit
  lo   <- uniroot(obj, c(est - 4 * se, est - se / 100))$root
  hi   <- uniroot(obj, c(est + se / 100, est + 4 * se))$root
  list(lower = lo, upper = hi, LR = lr, P = 1 - pchisq(lr, 1))
}

# 1. Right-censored single-covariate cph() across efron, breslow, exact ties
for(meth in c('efron', 'breslow', 'exact')) {
  f <- cph(Surv(futime, fustat) ~ rx, data = ovarian, method = meth,
           x = TRUE, y = TRUE)
  co <- contrast(f, list(rx = 2), list(rx = 1), conf.type = 'profile')
  ref <- exact_profile(Surv(futime, fustat) ~ rx, '1', ovarian, ovarian$rx,
                       ties = meth)
  stopifnot(
    abs(co$Lower  - ref$lower) < 1e-4,
    abs(co$Upper  - ref$upper) < 1e-4,
    abs(co$LR     - ref$LR)    < 1e-6,
    abs(co$Pvalue - ref$P)     < 1e-6
  )
}

# 2. Stratified single-covariate cph()
fs <- cph(Surv(futime, fustat) ~ rx + strat(ecog.ps), data = ovarian,
          x = TRUE, y = TRUE)
cos <- contrast(fs, list(rx = 2), list(rx = 1), conf.type = 'profile')
refs <- exact_profile(Surv(futime, fustat) ~ rx + strata(ecog.ps),
                      'strata(ecog.ps)', ovarian, ovarian$rx, ties = 'efron')
stopifnot(
  abs(cos$Lower  - refs$lower) < 1e-4,
  abs(cos$Upper  - refs$upper) < 1e-4,
  abs(cos$LR     - refs$LR)    < 1e-6,
  abs(cos$Pvalue - refs$P)     < 1e-6
)

# 3. Counting-process single-covariate cph() across efron, breslow, exact ties
ov_ag <- ovarian
ov_ag$tstart <- 0
dd_ag <- datadist(ov_ag)
options(datadist = 'dd_ag')
for(meth in c('efron', 'breslow', 'exact')) {
  fa <- cph(Surv(tstart, futime, fustat) ~ rx, data = ov_ag, method = meth,
            x = TRUE, y = TRUE)
  coa <- contrast(fa, list(rx = 2), list(rx = 1), conf.type = 'profile')
  refa <- exact_profile(Surv(tstart, futime, fustat) ~ rx, '1', ov_ag, ov_ag$rx,
                        ties = meth)
  stopifnot(
    abs(coa$Lower  - refa$lower) < 1e-4,
    abs(coa$Upper  - refa$upper) < 1e-4,
    abs(coa$LR     - refa$LR)    < 1e-6,
    abs(coa$Pvalue - refa$P)     < 1e-6
  )
}

# 4. Weighted single-covariate cph()
ov_w <- ovarian
ov_w$wt <- rep(c(1, 2), length.out = nrow(ov_w))
dd_w <- datadist(ov_w)
options(datadist = 'dd_w')
fw <- cph(Surv(futime, fustat) ~ rx, data = ov_w, weights = wt,
          x = TRUE, y = TRUE)
cow <- contrast(fw, list(rx = 2), list(rx = 1), conf.type = 'profile')
refw <- exact_profile(Surv(futime, fustat) ~ rx, '1', ov_w, ov_w$rx,
                      ties = 'efron', wvec = ov_w$wt)
stopifnot(
  abs(cow$Lower  - refw$lower) < 1e-4,
  abs(cow$Upper  - refw$upper) < 1e-4,
  abs(cow$LR     - refw$LR)    < 1e-6,
  abs(cow$Pvalue - refw$P)     < 1e-6
)

cat('All single-covariate cph profile likelihood tests passed.\n')
