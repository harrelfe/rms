# From depigner package https://github.com/CorradoLanera/depigner/blob/master/tests/testthat/test-summary_interact.R

require(rms)

data("transplant", package = "survival")
transplant <- transplant[transplant[["event"]] != "censored", ] |>
  droplevels()

dd <- datadist(transplant); options(datadist='dd')
# Note that event is being treated as ordinal, which may be a problem
# Default lrm with default tol=1e-13 finds singular matrix
# Adding either transx=TRUE or tol=1e-14 fixes this
f <- lrm(event ~ rcs(age, 3) * (sex + abo) + rcs(year, 3), data=transplant,
         transx=TRUE, trace=1, opt_method='NR')
f <- lrm(event ~ rcs(age, 3) * (sex + abo) + rcs(year, 3), data=transplant,
         tol=1e-14, trace=1, opt_method='NR')

