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

