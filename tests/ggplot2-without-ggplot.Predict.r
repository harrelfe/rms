require(rms)
require(ggplot2)
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
sex            <- factor(sample(c('female','male'), n, TRUE))
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
d <- data.frame(y, blood.pressure, age, cholesterol, sex)
rm(y, blood.pressure, age, cholesterol, sex)
dd <- datadist(d); options(datadist='dd')

f <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
         data=d)
p <- Predict(f, cholesterol)
class(p) <- setdiff(class(p), 'Predict')
a <- attributes(p)$info$Design
g <- ggplot(p, aes(x=cholesterol, y=yhat)) + geom_line()
xl <- labelPlotmath(a$label['blood.pressure'],
                    a$units['blood.pressure'])
xl2 <- labelPlotmath(a$label['cholesterol'],
                     a$units['cholesterol'])
g <- g + xlab(xl)
g <- g + geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2, linetype=0)
g
g + histSpikeg(yhat ~ cholesterol, p, d, ylim=c(-1, 1.25))
g + histSpikeg(yhat ~ cholesterol, data=d, ylim=c(-1, 1.25))
g + histSpikeg(yhat ~ cholesterol, data=d, ylim=c(-1, 1.25), side=3)

p <- Predict(f, cholesterol, sex)
class(p) <- setdiff(class(p), 'Predict')
g <- ggplot(p, aes(x=cholesterol, y=yhat, color=sex)) + geom_line() +
  xlab(xl2) + ylim(-1, 1)
# show.legend=FALSE gets rid of slash in legend boxes
# See http://stackoverflow.com/questions/10660775/ggplot-legend-slashes
g <- g + geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2,
                linetype=0, show.legend=FALSE)
g
g + histSpikeg(yhat ~ cholesterol + sex, p, d, ylim=c(-1, 1.25))

p <- Predict(f, sex)
class(p) <- setdiff(class(p), 'Predict')
ggplot(p, aes(x=sex, y=yhat)) + coord_flip() + geom_point() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=0)

p     <- Predict(f)
a     <- attributes(p)$info
yl    <- a$ylabPlotmath
xlabs <- a$Design$label
unts  <- a$Design$units
ylim <- range(pretty(
  if(TRUE) c(p$yhat, p$lower, p$upper)
  else p$yhat), na.rm=TRUE)

grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
nr <- 1; nc <- 0
for(w in unique(p$.predictor.)) {
  nc <- nc + 1
  if(nc > 2) {nr <- nr + 1; nc <- 1}
  i <- p$.predictor. == w
  z <- p[i, w]
  yhat <- p[i, 'yhat']
  l  <- levels(z)
  ll <- length(l)
  xl <- labelPlotmath(xlabs[w], unts[w])
  zz <- data.frame(z, yhat)
  g  <- ggplot(zz, aes(x=z, y=yhat)) + ylim(ylim) +
    theme(plot.margin = unit(rep(.2, 4), 'cm'))

  g <- g + if(ll) geom_point() else geom_line()
  g <- g + xlab(xl) + ylab(yl)
  g <- g + if(ll)
    geom_errorbar(data=p[i,], aes(ymin=lower, ymax=upper), width=0)
    else
      geom_ribbon(data=p[i,], aes(ymin=lower, ymax=upper), alpha=0.2,
                  linetype=0, show.legend=FALSE)
  print(g, vp = grid::viewport(layout.pos.row = nr, layout.pos.col = nc))
}

# Change y scale to be uniform
# try to narrow gaps



p <- Predict(f, age, sex, blood.pressure=c(120,140,160),
             cholesterol=c(180,200,215))
class(p) <- setdiff(class(p), 'Predict')
g <- ggplot(p, aes(x=age, y=yhat, color=sex)) + geom_line()
g <- g + geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2,
                linetype=0, show.legend=FALSE)
g + facet_grid(blood.pressure ~ cholesterol)
g + facet_grid(cholesterol ~ blood.pressure)
eval(parse(text='g + facet_grid(cholesterol ~ blood.pressure)'))


# attr(p, 'info')$varying shows 4 predictors varying in order: age bp ch sex

g <- ggplot(p, aes(x=age, y=yhat)) + geom_line()
g <- g + geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2,
                linetype=0, show.legend=FALSE)
g + facet_grid(blood.pressure ~ cholesterol*sex)
g + facet_grid(cholesterol*sex ~ blood.pressure)

# Add superposition
g <- ggplot(p, aes(x=age, y=yhat, color=blood.pressure)) + geom_line()
g <- g +
  geom_ribbon(data=p, aes(ymin=lower, ymax=upper), alpha=0.2,
                linetype=0, show.legend=FALSE)
g=g + facet_grid(sex ~ blood.pressure)

if(FALSE) {   # doesn't work - where is .predictor.?
p <- as.data.frame(p)
g <- ggplot(p, aes(y=yhat)) + facet_wrap(~ .predictor., scales='free_x') +
  xlab(NULL)
require(plyr)
pa <- subset(p, .predictor. == 'age')
pc <- subset(p, .predictor. == 'cholesterol')

g <- g + geom_line(subset=.(.predictor.=='age'), aes(x=age)) +
  geom_ribbon(subset=.(.predictor.=='age'), aes(x=age, ymin=lower, ymax=upper),
              alpha=0.2, linetype=0, show.legend=FALSE) +
  geom_line(subset=.(.predictor.=='cholesterol'), aes(x=cholesterol)) +
  geom_ribbon(subset=.(.predictor.=='cholesterol'),
              aes(x=cholesterol, ymin=lower, ymax=upper),
              alpha=0.2, linetype=0, show.legend=FALSE)
g
g + geom_point(subset=.(.predictor.=='sex'), aes(x=as.numeric(sex))) +
    geom_errorbar(subset=.(.predictor.=='sex'),
               aes(x=as.numeric(sex), ymin=lower, ymax=upper), width=0)


## Will not work:
##  g + geom_point(subset=.(.predictor.=='sex'), aes(x=sex)) +
##  geom_errorbar(subset=.(.predictor.=='sex'),
##                aes(x=sex, ymin=lower, ymax=upper), width=0)
## Error: Discrete value supplied to continuous scale

xx <- NULL
pred <- p$.predictor.
for(i in unique(pred)) xx <- c(xx, p[pred == i, i])
p$xx <- xx
z <- ggplot(p, aes(x=xx, y=yhat)) +
     facet_wrap(~ .predictor., scales='free_x') + xlab(NULL) +
     geom_line() + geom_ribbon(aes(x=xx, ymin=lower, ymax=upper),
                            alpha=0.2, linetype=0, show.legend=FALSE)
z
}


## From http://stackoverflow.com/questions/11979017/changing-facet-label-to-math-formula-in-ggplot2
facet_wrap_labeller <- function(gg.plot, labels=NULL) {
  require(gridExtra)
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs      
  strips <- grep("strip_t", names(gg))
  for(ii in seq_along(labels))  {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text", 
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <-
      editGrob(modgrob,label=labels[ii])
  }
  g$grobs <- gg
  class(g) = c("arrange", "ggplot", class(g)) 
  g
}
if(FALSE) {
pold <- p
p$.predictor. <- factor(p$.predictor., names(a$label))
pmlabels <- vector('expression', length(a$label))
names(pmlabels) <- levels(p$.predictor.)
for(v in names(a$label)) pmlabels[v] <- 
 labelPlotmath(a$label[v], a$units[v])

## Re-order panels by original model specification
z <- ggplot(p, aes(x=xx, y=yhat)) +
     facet_wrap(~ .predictor., scales='free_x', ncol=3) + xlab(NULL) +
     geom_line() + geom_ribbon(aes(x=xx, ymin=lower, ymax=upper),
                            alpha=0.2, linetype=0, show.legend=FALSE)

facet_wrap_labeller(z, pmlabels)
}
