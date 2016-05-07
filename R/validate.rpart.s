validate.rpart <- function(fit, method, B, bw, rule, type, sls, aics,
                           force, estimates, pr=TRUE, k, rand, xval = 10, 
                           FUN, ...)
{
  if(missing(FUN)) FUN <- function(..., k) rpart::prune(..., cp=k)

  act <- (fit$call)$na.action
  if(! length(act))  act <- function(x) x
  m <- model.frame(fit, na.action = act)
  if(! is.data.frame(m)) stop('you must specify model=T in the fit')
  y <- model.extract(m, 'response')
  ytype <- if(inherits(y, 'Surv')) 'Surv'
  else
    if(is.logical(y) ||
       ((length(un <- sort(unique(y[! is.na(y)]))) == 2) && un[1] == 0 &&
        un[2] == 1)) 'binary'
  else 'other'
  if(ytype == 'binary' && is.factor(y)) y <- as.numeric(y) - 1
  dxyf <- switch(ytype,
                 binary = function(x, y) somers2(x, y)['Dxy'],
                 Surv   = function(x, y) - dxy.cens(x, y)['Dxy'],
                 other  = function(x, y) dxy.cens(x, y)['Dxy'])

  call <- match.call()
  method <- call$method
  size <- NULL
  if(missing(k)) {
    k    <- fit$cptable[, 'CP']
    size <- fit$cptable[, 'nsplit']
  }
  if(missing(rand))
    rand <- sample(xval, NROW(m[[1]]), replace = TRUE)
  which <- unique(rand)
  pdyx.app <- pdyx.val <- pb.app <- pb.val <- double(length(k))
  l <- 0
  for(kk in k) {
    l <- l + 1
    dyx.val <- dyx.app <- b.val <- b.app <- double(length(which))
    j <- 0
    for(i in which) {
      j <- j + 1
      s <- rand != i
      tlearn <- rpart::rpart(model=m[s, ])
      papp <- if(kk == 0) tlearn else FUN(tlearn, k = kk, ...)
      if(nrow(papp$frame) == 1) {
        dyx.app[j] <- dyx.val[j] <- 0	#no splits
        if(ytype != 'Surv')
          b.app[j] <- b.val[j] <- mean((y - mean(y))^2, na.rm = TRUE)
      } else {
        yhat <- predict(papp, newdata = m[s,  ])
        if(is.matrix(yhat) && ncol(yhat) > 1)
          yhat <- yhat[, ncol(yhat), drop=TRUE]
        ysub <- if(ytype == 'Surv') y[s, ] else y[s]
        ## tree with factor binary y
        if(ytype != 'Surv') b.app[j] <- mean((yhat - ysub)^2)
        dyx.app[j] <- dxyf(yhat, ysub)
        s <- rand == i
        yhat <- predict(papp, newdata = m[s,  ])
        ysub <- if(ytype == 'Surv') y[s, ] else y[s]
        if(ytype != 'Surv') b.val[j] <- mean((yhat - ysub)^2)
        dyx.val[j] <- dxyf(yhat, ysub)
      }
    }
    pdyx.app[l] <- mean(dyx.app)
    pdyx.val[l] <- mean(dyx.val)
    pb.app[l]   <- mean(b.app)
    pb.val[l]   <- mean(b.val)
    if(pr) {
      dyx.app <- c(dyx.app, pdyx.app[l])
      dyx.val <- c(dyx.val, pdyx.val[l])
      b.app <- c(b.app, pb.app[l])
      b.val <- c(b.val, pb.val[l])
      cat("\n\nk=", format(kk), ":\n\n")
      rnam <- c(as.character(1 : j), "Mean")
      if(ytype == 'Surv') {
        dyx <- cbind(dyx.app, dyx.val)
        dimnames(dyx) <- list(rnam, c('Dxy Training', 'Dxy Test'))
      } else {
        dyx <- cbind(dyx.app, dyx.val, b.app, b.val)
        dimnames(dyx) <- list(rnam,
                              c("Dxy Training", "Dxy Test", "MSE Training", 
                                "MSE Test"))
      }
      print(dyx)
    }
  }
  if(ytype == 'Surv') pb.app <- pb.val <- NULL
  structure(list(k = k, size = size, dxy.app = pdyx.app,
                 dxy.val = pdyx.val, mse.app = pb.app, mse.val = pb.val,
                 ytype = ytype, xval = xval),
            class = "validate.rpart")
}

print.validate.rpart <- function(x, ...)
{
  cat(x$xval, "-fold cross-validation\n\n", sep = "")
  w <- cbind(k = x$k, size = x$size, Dxy.apparent = x$dxy.app, 
             Dxy.val = x$dxy.val, MSE.apparent = x$mse.app,
             MSE.val = x$mse.val)
  if(x$ytype == 'binary')
    dimnames(w) <- list(NULL, c("k", if(length(x$size)) "size", 
                                "Dxy.apparent", "Dxy.val", "Brier.apparent", 
                                "Brier.val"))
  invisible(print(w))
}

plot.validate.rpart <- function(x, what = c("mse", "dxy"),
                                legendloc = locator, ...)
{
  if(! missing(what) && x$ytype == 'Surv' && 'mse' %in% what)
    stop('may not specify what="dxy" for survival trees')
  if(x$ytype == 'Surv') what <- 'dxy'
  
  obj <- x
  if(length(obj$size)) {
    x <- obj$size
    xlab <- "Number of Nodes"
  } else {
    x <- obj$k
    xlab <- "Cost/Complexity Parameter"
  }
  if("mse" %in% what) {
    blab <- if(obj$ytype == 'binary') "Brier Score" else "Mean Squared Error"
    ylim <- range(c(obj$mse.app, obj$mse.val))
    plot(x, obj$mse.app, xlab = xlab, ylab = blab, ylim = ylim, 
         type = "n")
    lines(x, obj$mse.app, lty = 3)
    lines(x, obj$mse.val, lty = 1)
    title(sub = paste(obj$xval, "-fold cross-validation", sep = ""),
          adj = 0)
    if(is.function(legendloc))
      legend(legendloc(1), c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n")
    else {
      legend(grconvertX(legendloc[1], from='npc'),
             grconvertY(legendloc[2], from='npc'),
             c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n")
    }
  }
  if("dxy" %in% what) {
    ylim <- range(c(obj$dxy.app, obj$dxy.val))
    plot(x, obj$dxy.app, xlab = xlab, ylab = "Somers' Dxy",
         ylim = ylim, type = "n")
    lines(x, obj$dxy.app, lty = 3)
    lines(x, obj$dxy.val, lty = 1)
    title(sub = paste(obj$xval, "-fold cross-validation", sep = ""),
          adj = 0)
    if(is.function(legendloc))
      legend(legendloc(1), c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n")
    else {
      par(usr=c(0,1,0,1))
      legend(legendloc[1],legendloc[2],
             c("Apparent", "Cross-validated"),
             lty = c(3, 1), bty = "n")
    }
  }
  invisible()
}
