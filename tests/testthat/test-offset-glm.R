library(testthat)
library(rms)
# sessionInfo("rms") # v. 4.2-1 from github
test_df <- data.frame(Y  = c(15,  7, 36,  4, 16, 12, 41, 15),
                      N  = c(4949, 3534, 12210, 344, 6178, 
                             4883, 11256, 7125),
                      x1 = c(-0.1, 0, 0.2, 0, 1, 1.1, 1.1, 1),
                      x2 = c(2.2, 1.5, 4.5, 7.2, 4.5, 3.2, 9.1, 5.2))
test_df$fact_var <- factor(rep(c("A", "B", "C"), times=c(4,2, nrow(test_df)-6)))
ddist <- datadist(test_df)
options(datadist="ddist")

test_that("Compare basic Glm with glm",{
  # Base check
  fit.Glm <- Glm(Y ~ x1 + x2 + fact_var, data=test_df)
  fit.glm <- glm(Y ~ x1 + x2 + fact_var, data=test_df)
  expect_equivalent(coef(fit.Glm), coef(fit.glm))
  nd <- data.frame(x1=1, x2=4, fact_var="A")
  expect_equivalent(predict(fit.Glm, newdata=nd),
                    predict(fit.glm, newdata=nd))
  nd <- data.frame(x1=1, x2=4, fact_var="B")
  expect_equivalent(predict(fit.Glm, newdata=nd),
                    predict(fit.glm, newdata=nd))
  
  expect_equivalent(vcov(fit.Glm),
                    vcov(fit.glm),
                    info="The covariance matrices should be identical")
})

test_that("Compare offset Glm with glm with data parameter",{
  # Check offset term
  fit.Glm <- list(plain=list(formula = Glm(Y ~ x1 + x2 + fact_var + offset(N), 
                                           data=test_df),
                             parameter = Glm(Y ~ x1 + x2 + fact_var , offset=N, 
                                             data=test_df)),
                  poisson=list(formula = Glm(Y ~ x1 + x2 + fact_var + offset(log(N)), 
                                             data=test_df, family=poisson),
                               parameter = Glm(Y ~ x1 + x2 + fact_var , offset=log(N), 
                                               data=test_df, family=poisson)))
  fit.glm <- list(plain=glm(Y ~ x1 + x2 + fact_var , offset=N, 
                            data=test_df),
                  poisson=glm(Y ~ x1 + x2 + fact_var , offset=log(N), 
                              data=test_df, family=poisson))
  # Check offset term
  for (glm_family in names(fit.Glm)){
    for (fn in names(fit.Glm[[glm_family]])){
      expect_equivalent(coef(fit.Glm[[glm_family]][[fn]]), 
                        coef(fit.glm[[glm_family]]),
                        info=paste("Coefficients are not right",
                                   "for", fn,
                                   "for family =", glm_family))
      
      expect_equivalent(contrast(fit.Glm[[glm_family]][[fn]], 
                                 a=list(x1=1), 
                                 b=list(x1=0))$Contrast,
                        coef(fit.Glm[[glm_family]][[fn]])["x1"],
                        info=paste("Contrast fails to provide the correct coefficient",
                                   "for", fn,
                                   "for family =", glm_family))
      
      expect_equivalent(vcov(fit.Glm[[glm_family]][[fn]]),
                        vcov(fit.glm[[glm_family]]),
                        info=paste("The covariance matrices should be identical", 
                                   "when using", fn,
                                   "for family =", glm_family))
      
      
      nd <- data.frame(x1=1, x2=4, N=0, fact_var="A")
      expect_equivalent(predict(fit.Glm[[glm_family]][[fn]], newdata=nd),
                        predict(fit.glm[[glm_family]], newdata=nd),
                        info=paste("Offset term is not equal to 0",
                                   "for the", fn, "version",
                                   "for family =", glm_family))
      nd$N <- 1000
      expect_equivalent(predict(fit.Glm[[glm_family]][[fn]], newdata=nd),
                        predict(fit.glm[[glm_family]], newdata=nd),
                        info=paste("Offset term is not equal to 1000",
                                   "for the", fn, "version",
                                   "for family =", glm_family))
    }
  }
})

test_that("Compare offset Glm with lm with(dataset, Glm(...))",{
  fit.Glm <- list(plain=list(param = with(test_df, 
                                          Glm(Y ~ x1 + x2 + fact_var, offset=N)),
                             formula = with(test_df, 
                                            Glm(Y ~ x1 + x2 + offset(N) + fact_var))),
                  poisson=list(param = with(test_df, 
                                             Glm(Y ~ x1 + x2 + fact_var, offset=log(N), family=poisson)),
                                formula = with(test_df, 
                                               Glm(Y ~ x1 + x2 + offset(log(N)) + fact_var,, family=poisson))))
  fit.glm <- list(plain=with(test_df, 
                             glm(Y ~ x1 + x2 + fact_var, offset=N)),
                  poisson=with(test_df, 
                               glm(Y ~ x1 + x2 + fact_var, offset=log(N), family=poisson)))
                  
                  
  # Check offset term
  for (glm_family in names(fit.Glm)){
    for (fn in names(fit.Glm[[glm_family]])){
      expect_equivalent(coef(fit.Glm[[glm_family]][[fn]]), 
                        coef(fit.glm[[glm_family]]),
                        info=paste("Coefficients are not right",
                                   "for", fn,
                                   "for family =", glm_family))
      
      expect_equivalent(contrast(fit.Glm[[glm_family]][[fn]], 
                                 a=list(x1=1), 
                                 b=list(x1=0))$Contrast,
                        coef(fit.Glm[[glm_family]][[fn]])["x1"],
                        info=paste("Contrast fails to provide the correct coefficient",
                                   "for", fn,
                                   "for family =", glm_family))
      
      expect_equivalent(vcov(fit.Glm[[glm_family]][[fn]]),
                        vcov(fit.glm[[glm_family]]),
                        info=paste("The covariance matrices should be identical", 
                                   "when using", fn,
                                   "for family =", glm_family))
      
      
      nd <- data.frame(x1=1, x2=4, N=0, fact_var="A")
      expect_equivalent(predict(fit.Glm[[glm_family]][[fn]], newdata=nd),
                        predict(fit.glm[[glm_family]], newdata=nd),
                        info=paste("Offset term is not equal to 0",
                                   "for the", fn, "version",
                                   "for family =", glm_family))
      nd$N <- 1000
      expect_equivalent(predict(fit.Glm[[glm_family]][[fn]], newdata=nd),
                        predict(fit.glm[[glm_family]], newdata=nd),
                        info=paste("Offset term is not equal to 1000",
                                   "for the", fn, "version",
                                   "for family =", glm_family))
    }
  }
})