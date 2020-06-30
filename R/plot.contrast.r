##' Plot Bayesian Contrast Posterior Densities
##'
##' If there are exactly two contrasts and \code{bivar=TRUE} plots an elliptical or kernal (based on \code{bivarmethod} posterior density contour with probability \code{prob}).  Otherwise plots a series of posterior densities of contrasts along with HPD intervals, posterior means, and medians.  When the result being plotted comes from \code{contrast} with \code{fun=} specified, both the two individual estimates and their difference are plotted.
##' @title plot.contrast.rms
##' @param x the result of \code{contrast.rms}
##' @param bivar set to \code{TRUE} to plot 2-d posterior density contour
##' @param bivarmethod see \code{pdensityCountour}
##' @param prob posterior coverage probability for HPD interval or 2-d contour
##' @param which applies when plotting the result of \code{contrast(..., fun=)}, defaulting to showing the posterior density of both estimates plus their difference.  Set to \code{"ind"} to only show the two individual densities or \code{"diff"} to only show the posterior density for the differences.
##' @param nrow for \code{ggplot2::facet_wrap}
##' @param ncol likewise
##' 
##' @param ... unused
##' @return \code{ggplot2} object
##' @author Frank Harrell
plot.contrast.rms <- function(x, bivar=FALSE,
                              bivarmethod=c('ellipse', 'kernel'), prob=0.95,
                              which=c('both', 'diff', 'ind'),
                              nrow=NULL, ncol=NULL, ...) {
  bivarmethod <- match.arg(bivarmethod)
  which       <- match.arg(which)

  if('esta' %in% names(x)) {
    # Handle output pertaining to fun= result differently
    esta <- x$esta
    estb <- x$estb
    w <- function(draws, what) {
      nd     <- nrow(draws)
      theta  <- as.vector(draws)
      contr  <- colnames(draws)
      if(length(contr) == 1 && contr == '1') contr <- ''
      cont   <- factor(rep(contr, each=nd), contr)
      d      <- data.frame(what, contr=cont, theta)
      f <- function(x) {
        hpd <- HPDint(x, prob)
        r <- c(mean(x), median(x), hpd)
        names(r) <- c('Mean', 'Median', 'Lower', 'Upper')
        r
      }
      est  <- apply(draws, 2, f)
      stat <- rownames(est)
      stat   <- ifelse(stat %in% c('Lower', 'Upper'),
                       paste(prob, 'HPDI'), stat)
      eparam <- factor(rep(contr, each=nrow(est)), contr)
      stat   <- rep(stat, length(contr))
      est    <- as.vector(est)
      de <- data.frame(what, contr=eparam, est, stat)
      list(d=d, de=de)
    }

    w1 <- w(esta, 'First')
    w2 <- w(estb, 'Second')
    w3 <- w(esta - estb, 'First - Second')
    d  <- rbind(w1$d,  w2$d,  w3$d)
    de <- rbind(w1$de, w2$de, w3$de)
    lev <- c('First', 'Second', 'First - Second')
    d$what  <- factor(d$what,  lev)
    de$what <- factor(de$what, lev)
    if(which == 'diff') {
      d  <- subset(d, what == 'First - Second')
      de <- subset(de, what == 'First - Second')
    } else if(which == 'ind') {
      d  <- subset(d, what != 'First - Second')
      de <- subset(de, what != 'First - Second')
      }
      
    g <- ggplot(d, aes(x=theta)) + geom_density() +
      geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4))) +
      facet_grid(what ~ contr) +
      guides(color=guide_legend(title='')) +
      xlab('') + ylab('')
    return(g)
  }

  cdraws <- x$cdraws
  if(! length(cdraws))
    stop('plot method for contrast.rms objects implemented only for Bayesian models')
  nd <- nrow(cdraws)
  cn <- colnames(cdraws)
  if(all(cn == as.character(1 : ncol(cdraws))))
    cn <- paste('Contrast', cn)
  colnames(cdraws) <- cn
  
  if(ncol(cdraws) == 2 && bivar) {
    g <- pdensityContour(cdraws[, 1], cdraws[, 2], prob=prob, pl=TRUE,
                         method=bivarmethod)
    g <- g + xlab(cn[1]) +  ylab(cn[2])
    return(g)
  }

  hpd   <- apply(cdraws, 2, HPDint, prob=prob)
  
  draws  <- as.vector(cdraws)
  which  <- colnames(cdraws)
  param  <- factor(rep(which, each=nd), which)

  g      <- function(x) c(mean=mean(x), median=median(x))
  est    <- apply(cdraws, 2, g)
  est    <- rbind(est, hpd)

  stat   <- rownames(est)
  stat   <- ifelse(stat %in% c('Lower', 'Upper'),
                   paste(prob, 'HPDI'), stat)
  eparam <- factor(rep(which, each=nrow(est)), which)
  stat   <- rep(stat, length(which))
  est    <- as.vector(est)
 
  d  <- data.frame(param, draws)
  de <- data.frame(param=eparam, est, stat)
  g <- ggplot(d, aes(x=draws)) + geom_density() +
         geom_vline(data=de, aes(xintercept=est, color=stat, alpha=I(0.4))) +
         facet_wrap(~ param, scales='free', nrow=nrow, ncol=ncol) +
         guides(color=guide_legend(title='')) +
         xlab('') + ylab('')
  g
}

utils::globalVariables(c('theta', 'what'))
