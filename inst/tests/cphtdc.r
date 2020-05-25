# From Max Gordon
require(rms)
# Same simulated data set as used in the cph-example
n <- 1000
set.seed(731)
age <- 50 + 12*rnorm(n)
sex <- factor(sample(c('Male','Female'), n, 
                     rep=TRUE, prob=c(.6, .4)))
cens <- 15*runif(n)
h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
dt <- -log(runif(n))/h
e <- ifelse(dt <= cens,1,0)
dt <- pmin(dt, cens)

test <- data.frame(age = age,
                      sex = sex,
                      Start = 0,
                      dt = dt,
                      e = e)

dd <<- datadist(test); options(datadist='dd')

f <- cph(Surv(dt,e) ~ rcs(age,4) + sex, x=TRUE, y=TRUE)
cox.zph(f, "rank")             # tests of PH
 

# Now to the actual time-interaction
if(! require(Epi)) quit(save='no')

lxs <- Lexis(entry = list(Timeband = Start),
             exit = list(Timeband = dt, Age = age + dt),
             exit.status = e,
             data = test)
subset(lxs, lex.id %in% 1:3)
  
spl <- 
  splitLexis(lxs, 
             time.scale = "Timeband",
             breaks = seq(from = 0, 
                          to = ceiling(max(lxs$lex.dur)),
                          by = .5))

subset(spl, lex.id %in% 1:3)

spl$Stop <- spl$Timeband + spl$lex.dur

dd <- datadist(spl)

#######################
# Regular interaction #
#######################
  coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex*Timeband,
      data = spl)
# Gives:
# Call:
# coxph(formula = Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
#         Age + sex * Timeband, data = spl)
# 
# 
# coef exp(coef) se(coef)     z       p
# Age               0.0420     1.043  0.00558  7.53 5.0e-14
# sexMale          -0.9457     0.388  0.26014 -3.64 2.8e-04
# Timeband              NA        NA  0.00000    NA      NA
# sexMale:Timeband  0.0868     1.091  0.05360  1.62 1.1e-01
# 
# Likelihood ratio test=72.7  on 3 df, p=1.11e-15  n= 13421, number of events= 183 
# Warning message:
#   In coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~  :
#              X matrix deemed to be singular; variable 3  

cph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex*Timeband,
    data = spl)
# Gives:
# X matrix deemed to be singular; variable Timeband 
# 
# Model Did Not Converge.  No summary provided.

###############################
# Forced singular interaction #
###############################
coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
        Age + sex + sex:Timeband,
      data = spl)
# Gives:
# Call:
# coxph(formula = Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
#         Age + sex + sex:Timeband, data = spl)
# 
# 
# coef exp(coef) se(coef)     z       p
# Age                 0.0420     1.043  0.00558  7.53 5.0e-14
# sexMale            -0.9457     0.388  0.26014 -3.64 2.8e-04
# sexFemale:Timeband -0.0868     0.917  0.05360 -1.62 1.1e-01
# sexMale:Timeband        NA        NA  0.00000    NA      NA
# 
# Likelihood ratio test=72.7  on 3 df, p=1.11e-15  n= 13421, number of events= 183 
# Warning message:
#   In coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~  :
#              X matrix deemed to be singular; variable 4

coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex + I((sex == "Male")*Timeband),
    data = spl)
# Gives:
# Call:
# coxph(formula = Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
#           Age + sex + I((sex == "Male") * Timeband), data = spl)
# 
# 
# coef exp(coef) se(coef)     z       p
# Age                            0.0420     1.043  0.00558  7.53 5.0e-14
# sexMale                       -0.9457     0.388  0.26014 -3.64 2.8e-04
# I((sex == "Male") * Timeband)  0.0868     1.091  0.05360  1.62 1.1e-01

cph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex + sex:Timeband,
    data = spl)
# Gives:
# X matrix deemed to be singular; variable sex=Male * NA 
# 
# Model Did Not Converge.  No summary provided. 

cph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
      rcs(Age, 4) + sex + asis((sex == "Male")*Timeband),
    data = spl)
# Gives:
# Err. in limits[[zname]] <- if (any(Limnames == zname)) { : 
#   more elements supplied than there are to replace

#############
# After fix #
#############
fit_coxph <- coxph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex + I((sex == "Male")*Timeband),
                        data = spl)

fit_cph <- cph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ Age + sex + asis((sex == "Male")*Timeband),
    data = spl)

# Basically the same
cbind(coxph=coef(fit_coxph), 
      cph=coef(fit_cph))
# Although numerically not equivalent
expect_true(sum(abs(coef(fit_cph) - coef(fit_coxph))) < 
              .Machine$double.eps)

#############
# Needs fix #
#############

Predict(fit_cph)
#  Err. in asis((sex == "Male") * Timeband) : object 'Timeband' not found 

# Not really working as expected
contrast(fit_cph, 
         a=list(sex = "Male"),
         b=list(sex = "Female"))
#  Err. in Getlimi(name[i], Limval, need.all = TRUE) : 
#    no limits defined by datadist for variable sex_Timeband 
contrast(fit_cph, 
         a=list(sex = "Male",
                Timeband = 0),
         b=list(sex = "Female",
                Timeband = seq(0, 10, by=.1)))
# Err. in gendata(list(coefficients = c(0.0420352254526414, -0.945650117874665,  : 
#   factor(s) not in design: Timeband 

#Ok, thank you. I can get around the problem by manually generating an interaction variable - seems to work satisfactory:

spl_alt <- 
  within(spl, {
    Male_time_int = (sex == "Male")*Timeband
  })

spl_alt$lex.Cst <- NULL
spl_alt$Start <- NULL
dd <- datadist(spl_alt)
options(datadist = "dd")

model <-
  cph(Surv(time = Timeband, time2 = Stop, event = lex.Xst) ~ 
        Age + sex + Male_time_int,
      data = spl_alt)

contrast(model, 
         a = list(sex = "Male", Male_time_int = 0:5),
         b = list(sex = "Female", Male_time_int = 0))
