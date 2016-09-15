# From Yong Hao Pua <puayonghao@gmail.com>
require(rms)
n.subjects <- 100    # original number is 30 on the rms documentation
ages <- rnorm(n.subjects, 50, 15)
sexes <- factor(sample(c('female','male'), n.subjects, TRUE))
logit <- (ages - 50) / 5
prob <- plogis(logit) # true prob not related to sex
id <- sample(1:n.subjects, 300, TRUE) # subjects sampled multiple times
length(unique(id))   # note that don't always get n.subjects sampled
length(table(id))
table(table(id)) # frequencies of number of obs/subject
age <- ages[id]
sex <- sexes[id]
# In truth, observations within subject are independent:
y <- ifelse(runif(300) <= prob[id], 1, 0)
f <- lrm(y ~ lsp(age, 50) * sex, x=TRUE, y=TRUE )
g <- robcov(f, id)
g$clusterInfo

