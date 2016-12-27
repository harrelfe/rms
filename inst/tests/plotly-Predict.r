require(rms)
require(plotly)
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
sex            <- factor(sample(c('female','male'), n, TRUE))
country        <- factor(sample(c('US', 'Canada'),  n, TRUE))
i <- sex == 'female'
cholesterol <- numeric(n)
cholesterol[i]   <- rnorm(sum(i),   170, 15)
cholesterol[! i] <- rnorm(sum(! i), 200, 25)
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'
# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
     (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)
cholesterol[1:3] <- NA   # 3 missings, at random
d <- data.frame(y, blood.pressure, age, cholesterol, sex, country)
rm(y, blood.pressure, age, cholesterol, sex, country)
dd <- datadist(d); options(datadist='dd')

f <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)) + country,
         data=d)


p <- Predict(f, cholesterol, sex)
source('~/R/Hmisc/R/histSpikeg.s')
source('~/R/rms/R/plotp.Predict.s')
source('~/R/Hmisc/R/scat1d.s')
# plotp(p, rdata=d, ylim=c(-1,2))

i <- attr(p, 'info')
cllab <- if(i$conf.int) paste0(i$conf.int, ' C.L.')

class(p) <- setdiff(class(p), 'Predict')

fm <- function(x) format(x, digits=4)

# pm <- subset(p, sex == 'male')
a <- i$Design
bpl <- labelPlotmath(a$label['blood.pressure'],
                     a$units['blood.pressure'], html=TRUE)
chl <- labelPlotmath(a$label['cholesterol'],
                     a$units['cholesterol'], html=TRUE)
agl <- labelPlotmath(a$label['age'],
                     a$units['age'], html=TRUE)

a <- plot_ly()
ht <- with(p, paste0('cholesterol=', fm(cholesterol), '<br>', 
                       fm(yhat), ' [', fm(lower), ',', fm(upper), ']'))
j <- which(p$cholesterol == min(p$cholesterol))
ht[j] <- paste0(ht[j], '<br>Adjusted to:<br>', i$adjust[1])

a <- add_lines(a, data=p, x=~cholesterol, y=~yhat, color=~sex,
               text=~ht, hoverinfo='text')
a <- add_ribbons(a, data=p, x=~cholesterol, ymin=~lower, ymax=~upper,
                 color=~sex, hoverinfo='none')
source('~/R/Hmisc/R/histSpikeg.s')
a <- histSpikeg(yhat ~ cholesterol + sex, predictions=p,
                data=d, plotly=a, ylim=c(-1, 2))
layout(a, xaxis=list(title=chl),
          yaxis=list(title=i$ylabhtml, range=c(-1, 2)))

p <- Predict(f)
# w <- plotp(p, rdata=d)
# w$Continuous
# w$Categorical

i <- attr(p, 'info')
ylim <- range(c(p$lower, p$upper, p$yhat), na.rm=TRUE)
p <- subset(p, .predictor. %nin% c('sex', 'country'))
class(p) <- 'data.frame'
r <- subset(p, .predictor. == 'age')
r$ht <- with(r, paste0('age=', fm(age), '<br>', 
                       fm(yhat), ' [', fm(lower), ',', fm(upper), ']'))
r$ht[1] <- paste0(r$ht[1], '<br>Adjusted to:<br>', i$adjust[3])
a <- plot_ly(r)
a <- add_lines(a, x=~age, y=~yhat, text=~ht, color=I('black'), hoverinfo='text',
               name='yhat', legendgroup='yhat')
a <- add_ribbons(a, x=~age, ymin=~lower, ymax=~upper, color=I('lightgray'),
                 hoverinfo='none', name=cllab, legendgroup=cllab)
source('~/R/Hmisc/R/histSpikeg.s')
a <- histSpikeg(yhat ~ age, data=d, predictions=r, ylim=ylim, plotly=a)
#aa <- histSpikep(a, x=d$age, y=approx(r$age, r$yhat, xout=d$age)$y, z=1)
ex <- function(x, delta=0) {
    r <- range(x, na.rm=TRUE)
    if(delta == 0) return(r)
    c(r[1] - delta * diff(r), r[2] + delta * diff(r))
}
a <- plotly::layout(a, xaxis=list(title=agl, range=ex(d$age)))

r <- subset(p, .predictor. == 'cholesterol')
r$ht <- with(r, paste0('cholesterol=', fm(cholesterol), '<br>',
                       fm(yhat), ' [', fm(lower), ',', fm(upper), ']'))
r$ht[1] <- paste0(r$ht[1], '<br>Adjusted to:<br>', i$adjust[4])
b <- plot_ly(r)
b <- add_lines(b, x=~cholesterol, y=~yhat, text=~ht, color=I('black'), hoverinfo='text',
               name='yhat', showlegend=FALSE, legendgroup='yhat')
b <- add_ribbons(b, x=~cholesterol, ymin=~lower, ymax=~upper, color=I('lightgray'),
                 hoverinfo='none', name=cllab, showlegend=FALSE, legendgroup=cllab)
b <- histSpikeg(yhat ~ cholesterol, data=d, predictions=r, ylim=ylim, 
                plotly=b, showlegend=FALSE)
b <- layout(b, xaxis=list(title='cholesterol', range=ex(d$cholesterol)))
plotly::subplot(a, b, nrows=1, shareY=TRUE, titleX=TRUE)

p <- Predict(f)
r <- subset(p, .predictor. == 'sex')
a <- plot_ly(r, color=I('black'), height=plotlyParm$heightDotchart(2))
a <- add_segments(a, y=~sex, x=~lower, yend=~sex, xend=~upper, color=I('lightgray'), name=cllab, legendgroup=cllab)
a <- add_markers(a, y=~sex, x=~yhat, name='Estimate', legendgroup='Estimate')
#lm <- plotlyParm$lrmargin('female')
a <- layout(a, xaxis=list(title=i$ylabhtml), 
               yaxis=list(title='Sex', titlefont=list(size=10)))

r <- subset(p, .predictor. == 'country')
b <- plot_ly(r, color=I('black'), height=plotlyParm$heightDotchart(2))
b <- add_segments(b, y=~country, x=~lower, yend=~country, xend=~upper, color=I('lightgray'), name=cllab, legendgroup=cllab, showlegend=FALSE)
b <- add_markers(b, y=~country, x=~yhat, name='Estimate', legendgroup='Estimate', showlegend=FALSE)
#lm <- plotlyParm$lrmargin('Canada')
b <- layout(b, xaxis=list(title=i$ylabhtml),
               yaxis=list(title='Country', titlefont=list(size=10)))

plotly::subplot(a, b, shareX=TRUE, titleY=TRUE, nrows=2, heights=c(2, 2) / sum(c(2, 2)))

p <- Predict(f, sex)
class(p) <- setdiff(class(p), 'Predict')


p <- Predict(f, age, sex, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215))
