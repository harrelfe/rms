require(rms)
for (i in unique(pbcseq$id)) {
  if (i == 1) {
  l <- length(pbcseq$id[pbcseq$id==i])
  start <- pbcseq$day[pbcseq$id==i]
  stop <- c(pbcseq$day[pbcseq$id==i][2:l],pbcseq$futime[pbcseq$id==i][1])
  event <- c(rep(0,l-1),pbcseq$status[pbcseq$id==i][1])
  } else  {
   l <- length(pbcseq$id[pbcseq$id==i])
   if (l==1) {
   t1 <- pbcseq$day[pbcseq$id==i]
   t2 <- pbcseq$futime[pbcseq$id==i]
   e <- pbcseq$status[pbcseq$id==i]
   start <- c(start,t1)
   stop <- c(stop,t2)
   event <- c(event,e)
   } else if (l>1) {
   t1 <- pbcseq$day[pbcseq$id==i]
   t2 <- c(pbcseq$day[pbcseq$id==i][2:l],pbcseq$futime[pbcseq$id==i][1])
   e <- c(rep(0,l-1),pbcseq$status[pbcseq$id==i][1])
   start <- c(start,t1)
   stop <- c(stop,t2)
   event <- c(event,e)
   }
 }
}

pbcseq <- data.frame(pbcseq,start,stop,event)

#bili is time-dependent covariate
fit  <-  cph(Surv(start, stop, event==2) ~ sex + log(bili) + rcs(age, 4),
             surv=T, x=T,y=T, data=pbcseq, eps=1e-8)

temp <- pbcseq[1:2,] #First id
temp$S <- with(temp, Surv(start, stop, event==2))
surv1 <- survfit(fit, newdata=temp, individual=TRUE)
surv2 <- survest(fit, newdata=temp, individual=TRUE)

A <- with(pbcseq, rcspline.eval(age, nk=4, inclx=TRUE))
temp$A <- A[1:2, ]

f  <-  coxph(Surv(start, stop, event==2) ~ sex + log(bili) + A,
             x=TRUE, y=TRUE, data=pbcseq)
s1 <- survfit(f, newdata=temp, individual=TRUE)
