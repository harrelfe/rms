# From Yong Hao Pua <puayonghao@gmail.com>
require(rms)
n<-20
set.seed(88)
age <- rnorm(n, 50, 10)
height <- rnorm(n, 1.7, 0.5)
cholesterol <- rnorm(n, 200, 25)
ch <- cut2(cholesterol, g=40, levels.mean=TRUE)
sex <- factor(sample(c("female","male"), n,TRUE))
dbase= data.frame(sex, age, height, cholesterol, ch)
dbase.dd <- datadist(dbase)
options(datadist = "dbase.dd")
fit <- ols (cholesterol ~ sex + height + age, x=T, y=T, data=dbase)
pentrace(fit, seq(0, 20, by = 0.1))
