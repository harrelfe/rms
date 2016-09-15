# From Pedro Emmanuel Alvarenga Americano do Brasil emmanuel.brasil@gmail.com

require(rms)
set.seed(1)
n  <- 400
n1 <- 300; n2 <- 100
data <- data.frame(outcome=c(rnorm(n1, mean = .052, sd = .005),
                     rnorm(n2, mean = .06, sd = .01)),
                   v2=sample(seq(20,80,5),n,T), 
                   v3=sample(seq(60,150,1),n,T), 
                   v4=c(rnorm(n1, mean = 80, sd = 10),
                     rnorm(n2, mean = 60, sd = 15)), 
                   v5=sample(c('M','M','F'),n,T), 
                   v6=c(rnorm(n1, mean = 80, sd = 10),
                     rnorm(n2, mean = 120, sd = 30))) 

# checking data
head(data)

# setting datadist
dd <- datadist(data);  options(datadist="dd")

# generating missings
m <- function() sample(1:n, 20, FALSE)

data$v2[m()] <- NA
data$v3[m()] <- NA
data$v4[m()] <- NA
data$v5[m()] <- NA
data$v6[m()] <- NA
plot(naclus(data))

# Imputing
imp <- aregImpute(~ outcome + v2 + v3 + v4 + v5 + v6, data, n.impute=10)

# fitting 
f <- fit.mult.impute(outcome ~ v6 + v2 + rcs(v3) + v5 * rcs(v4),
                     ols, imp, data, fit.reps=TRUE)
coef(f)
w <- NULL
for(i in 1 : 10) w <- rbind(w, coef(f$fits[[i]]))
w

s <- summary(f)
s
unclass(s)   # Effects are non-zero but small
plot(s)
