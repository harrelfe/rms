require(rms)

#Get a dataset/keep a few columns
load('boys.rda')   # originally from mice package
d <- boys[,c("age", "bmi", "reg")]
i <- with(d, is.na(bmi) | is.na(reg))
length(unique(d$age)); length(unique(d$age[! i]))

#sum(is.na(dat$age)) #0
####Models
##1) Complete case
#Set datadist
# dat_naomit <- na.omit(dat)
# dd <- datadist(dat_naomit)
# options(datadist = "dd")

dd <- datadist(d); options(datadist='dd')

#Run model
f <- orm(age ~ bmi + reg, data = d)

#Run a simple contrast
contrast(f, list(bmi = 20), list(bmi = 19))
summary(f, bmi=19:20, est.all=FALSE)

##2) Multiple imputation (default settings)

#Fit imputation model
# imp_mod <- mice(dat, m = 5) #Happens with ‘aregImpute’ as well
#Fit same orm model with imputed datasets
a <- aregImpute(~ age + bmi + reg, data=d, n.impute=5)

g <-
    fit.mult.impute(
        formula = age ~ bmi + reg,
        fitter = orm,
        xtrans = a,
        data = d
    )
dim(vcov(f, regcoef.only=TRUE))
dim(vcov(g, regcoef.only=TRUE))

summary(g, bmi=19:20, est.all=FALSE)

#Try the same contrast
contrast(g, list(bmi = 20), list(bmi = 19)) #Non-conformable dimension for matrix multiplication
