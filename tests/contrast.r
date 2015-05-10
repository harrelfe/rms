# From JoAnn Alvarez

require(rms)
set.seed(1)
age <- rnorm(200,40,12)
sex <- factor(sample(c('female','male'),200,TRUE))
logit <- (sex=='male') + (age-40)/5
y <- ifelse(runif(200) <= plogis(logit), 1, 0)
f <- lrm(y ~ pol(age,2)*sex)
anova(f)
# Compare a 30 year old female to a 40 year old male
# (with or without age x sex interaction in the model)
contrast(f, list(sex='female', age=30), list(sex='male', age=40))
# Test for interaction between age and sex, duplicating anova
contrast(f, list(sex='female', age=30),
            list(sex='male',   age=30),
            list(sex='female', age=c(40,50)),
            list(sex='male',   age=c(40,50)), type='joint')
# Duplicate overall sex effect in anova with 3 d.f.
contrast(f, list(sex='female', age=c(30,40,50)),
            list(sex='male',   age=c(30,40,50)), type='joint')


jim <- contrast(f, list(sex = "male", age=30),
                   list(sex = "male", age=40))
print(jim, fun = exp)
jane <- contrast(f, list(sex = c("male", "female"), age=30),
                    list(sex = c("male", "female"), age=40))
print(jane, fun = exp)




