# From Andy Bush <andy@kb4lsn.net>
require(rms)
set.seed(20)
x1<-10*runif(20,0,1)
y1<-c(rep(0,10),rep(1,10))
y2<-5*rnorm(20,0,1)

d<-data.frame(cbind(y1,y2,x1))
dd<-datadist(d)
options(datadist='dd')
flrm<-lrm(y1~x1,x=T,y=T,model=T)
nomlrm<-nomogram(flrm)
plot(nomlrm,xfac=.45)
fols<-ols(y2~x1,x=T,y=T,model=T)
nomols<-nomogram(fols)
plot(nomols,xfac=.45)


## From Zongheng Zhang zh_zhang1984@hotmail.com

n <- 1000    # sample size
set.seed(88) # set seed for replication
age<- rnorm(n, 65, 11)
lac<- round(abs(rnorm(n, 3, 1)),1)
sex<- factor(sample(1:2,n,prob=c(0.6,0.4),TRUE),
                      labels=c('male','female'))
shock<-factor(sample(1:4,n,prob=c(0.3,0.3,0.25,0.15),TRUE),
                labels=c('no','mild','moderate','severe'))
z<- 0.2*age + 3*lac* as.numeric(sex)+ 5*as.numeric(shock) -rnorm(n,36,15)
## linear combination with a bias

y <- ifelse(runif(n) <= plogis(z), 1, 0)
library(rms)
ddist <- datadist(age, lac, shock, sex)
options(datadist='ddist')
mod <- lrm(y ~ shock+lac*sex+age)
nom <- nomogram(mod,
        lp.at=seq(-3,4,by=0.5),
        fun=plogis,
        fun.at=c(.001,.01,.05,seq(.1,.9,by=.1),.95,.99,.999),
        funlabel="Risk of Death",
        conf.int=c(0.1, 0.7),
        abbrev=TRUE, #had not been working for shock
        minlength=1)

plot(nom, lplabel="Linear Predictor",
          fun.side=c(3,3,1,1,3,1,3,1,1,1,1,1,3),
          col.conf=c('red','green'),
          conf.space=c(0.1,0.5),
          label.every=3,
          col.grid = gray(c(0.8, 0.95)),
          which="shock")
legend.nomabbrev(nom, which='shock', x=.5, y=.5)
