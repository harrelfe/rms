library(testthat)
library(rms)
# sessionInfo("rms") # v. 4.2-1 from github
test_df <- data.frame(Y  = c(15,  7, 36,  4, 16, 12, 41, 15),
                      N  = c(4949, 3534, 12210, 344, 6178, 
                             4883, 11256, 7125),
                      x1 = c(-0.1, 0, 0.2, 0, 1, 1.1, 1.1, 1),
                      x2 = c(2.2, 1.5, 4.5, 7.2, 4.5, 3.2, 9.1, 5.2))
test_df$fact_var <- factor(rep(c("A", "B", "C"), times=c(4,2, nrow(test_df)-6)))
ddist <<- datadist(test_df)
options(datadist="ddist")

test_that("Compare basic ols with lm",{
  # Base check
  fit.ols <- ols(Y ~ x1 + x2 + fact_var, data=test_df)
  fit.lm <- lm(Y ~ x1 + x2 + fact_var, data=test_df)
  expect_equivalent(coef(fit.ols), coef(fit.lm))
  nd <- data.frame(x1=1, x2=4, fact_var="A")
  expect_equivalent(predict(fit.ols, newdata=nd),
                    predict(fit.lm, newdata=nd))
  nd <- data.frame(x1=1, x2=4, fact_var="B")
  expect_equivalent(predict(fit.ols, newdata=nd),
                    predict(fit.lm, newdata=nd))
  
  expect_equivalent(vcov(fit.ols),
                    vcov(fit.lm),
                    info="The covariance matrices should be identical")
})

test_that("Compare offset ols with lm with data parameter",{
  # Check offset term
  fit.ols <- list(formula = ols(Y ~ x1 + x2 + fact_var + offset(N), data=test_df),
                  parameter = ols(Y ~ x1 + x2 + fact_var , offset=N, data=test_df))
  fit.lm <- lm(Y ~ x1 + x2 + fact_var , offset=N, data=test_df)
  for (fn in names(fit.ols)){
    expect_equivalent(coef(fit.ols[[fn]]), 
                      coef(fit.lm),
                      info=paste("The coefficients are not the same",
                                 "when the ols fit is specified using",
                                 "an offset", fn))
    
    expect_equivalent(contrast(fit.ols[[fn]], 
                               a=list(x1=1), 
                               b=list(x1=0))$Contrast,
                      coef(fit.ols[[fn]])["x1"],
                      info=paste("Contrast fails to provide the correct coefficient",
                                 "for", fn))
    
    expect_equivalent(vcov(fit.ols[[fn]]),
                      vcov(fit.lm),
                      info=paste("The covariance matrices should be identical", 
                                 "when using", fn))
    
    
    nd <- data.frame(x1=1, x2=4, N=0, fact_var="A")
    expect_equivalent(predict(fit.ols[[fn]], newdata=nd),
                      predict(fit.lm, newdata=nd),
                      info=paste("The predictions are not the same",
                                 "for offset parameter=0 and",
                                 "when the ols fit is specified using",
                                 "an offset", fn))
    nd$N <- 1000
    expect_equivalent(predict(fit.ols[[fn]], newdata=nd),
                      predict(fit.lm, newdata=nd),
                      info=paste("The predictions are not the same",
                                 "for offset parameter=1000 and",
                                 "when the ols fit is specified using",
                                 "an offset", fn))
                      
  }
})

test_that("Compare offset ols with lm with(dataset, ols(...))",{
  fit.ols <- list(param = with(test_df, 
                               ols(Y ~ x1 + x2 + fact_var, offset=N)),
                  formula = with(test_df, 
                                 ols(Y ~ x1 + x2 + offset(N) + fact_var)))
  fit.lm <- with(test_df, 
                 lm(Y ~ x1 + x2 + fact_var, offset=N))

  # Check offset term
  for (fn in names(fit.ols)){
    expect_equivalent(coef(fit.ols[[fn]]), coef(fit.lm))

    expect_equivalent(contrast(fit.ols[[fn]], 
                               a=list(x1=1), 
                               b=list(x1=0))$Contrast,
                      coef(fit.ols[[fn]])["x1"],
                      info=paste("Contrast fails to provide the correct coefficient",
                                 "for", fn))
    
    expect_equivalent(vcov(fit.ols[[fn]]),
                      vcov(fit.lm),
                      info=paste("The covariance matrices should be identical", 
                                 "when using", fn))
    
    
    nd <- data.frame(x1=1, x2=4, N=0, fact_var="A")
    expect_equivalent(predict(fit.ols[[fn]], newdata=nd),
                      predict(fit.lm, newdata=nd),
                      info=paste("Offset term is not equal to 0",
                                 "for the", fn, "version"))
    nd$N <- 1000
    expect_equivalent(predict(fit.ols[[fn]], newdata=nd),
                      predict(fit.lm, newdata=nd),
                      info=paste("Offset term is not equal to 1000",
                                 "for the", fn, "version"))
  }
})

