validate.rpart <- function(fit, method, B, bw, rule, type, sls, aics,
                           force, estimates, pr=TRUE, k, rand, xval = 10, 
                           FUN, ...)
{
  if(missing(FUN))
    {
      require(rpart)
      FUN <- function(..., k) prune(..., cp=k)
    }
  act <- (fit$call)$na.action
  if(!length(act))
    act <- function(x) x
  m <- model.frame(fit, na.action = act)
  if(!is.data.frame(m)) stop('you must specify model=T in the fit')
  y <- model.extract(m, 'response')
  binary <- is.logical(y) ||
   ((length(un <- sort(unique(y[!is.na(y)]))) == 
    2) && un[1] == 0 && un[2] == 1)
  if(binary && is.factor(y)) y <- as.numeric(y) - 1
  call <- match.call()
  method <- call$method
  size <- NULL
  if(missing(k))
    {
      k    <- fit$cptable[,'CP']
      size <- fit$cptable[,'nsplit']
    }
  if(missing(rand))
    rand <- sample(xval, length(m[[1]]), replace = TRUE)
  which <- unique(rand)
  pdyx.app <- pdyx.val <- pb.app <- pb.val <- double(length(k))
  l <- 0
  for(kk in k)
    {
      l <- l + 1
      dyx.val <- dyx.app <- b.val <- b.app <- double(length(which))
      j <- 0
      for(i in which)
        {
          j <- j + 1
          s <- rand != i
          tlearn <- rpart(model=m[s,])
          papp <- if(kk == 0) tlearn else FUN(tlearn, k = kk, ...)
          if(nrow(papp$frame) == 1)
            {
              dyx.app[j] <- dyx.val[j] <- 0	#no splits
              b.app[j] <- b.val[j] <- mean((y - mean(y))^2, na.rm = TRUE)
            }
          else
            {
              yhat <- predict(papp, newdata = m[s,  ])
              if(is.matrix(yhat) && ncol(yhat) > 1)
                yhat <- yhat[,ncol(yhat),drop=TRUE]
              ## tree with factor binary y
              b.app[j] <- mean((yhat - y[s])^2)
              dyx.app[j] <- if(binary) somers2(yhat, y[s])["Dxy"] else
              rcorr.cens(yhat, y[s])["Dxy"]
              s <- rand == i
              yhat <- predict(papp, newdata = m[s,  ])
              b.val[j] <- mean((yhat - y[s])^2)
              dyx.val[j] <- if(binary) somers2(yhat, y[s])["Dxy"] else
              rcorr.cens(yhat, y[s])["Dxy"]
            }
        }
      pdyx.app[l] <- mean(dyx.app)
      pdyx.val[l] <- mean(dyx.val)
      pb.app[l]   <- mean(b.app)
      pb.val[l]   <- mean(b.val)
      if(pr)
        {
          dyx.app <- c(dyx.app, pdyx.app[l])
          dyx.val <- c(dyx.val, pdyx.val[l])
          b.app <- c(b.app, pb.app[l])
          b.val <- c(b.val, pb.val[l])
          cat("\n\nk=", format(kk), ":\n\n")
          dyx <- cbind(dyx.app, dyx.val, b.app, b.val)
          dimnames(dyx) <- list(c(as.character(1:j), "Mean"),
                                c("Dxy Training", "Dxy Test", "MSE Training", 
                                  "MSE Test"))
          print(dyx)
        }
    }
  structure(list(k = k, size = size, dxy.app = pdyx.app, dxy.val = 
                 pdyx.val, mse.app = pb.app, mse.val = pb.val,
                 binary = binary, xval = xval),
            class = "validate.rpart")
}

print.validate.rpart <- function(x, ...)
{
  cat(x$xval, "-fold cross-validation\n\n", sep = "")
  w <- cbind(k = x$k, size = x$size, Dxy.apparent = x$dxy.app, 
             Dxy.val = x$dxy.val, MSE.apparent = x$mse.app,
             MSE.val = x$mse.val)
  if(x$binary)
    dimnames(w) <- list(NULL, c("k", if(length(x$size)) "size", 
                                "Dxy.apparent", "Dxy.val", "Brier.apparent", 
                                "Brier.val"))
  invisible(print(w))
}

plot.validate.rpart <- function(x, what = c("mse", "dxy"),
                               legendloc = locator, ...)
{
  obj <- x
  if(length(obj$size))
    {
      x <- obj$size
      xlab <- "Number of Nodes"
    }
  else
    {
      x <- obj$k
      xlab <- "Cost/Complexity Parameter"
    }
  if("mse" %in% what)
    {
      blab <- if(obj$binary) "Brier Score" else "Mean Squared Error"
      ylim <- range(c(obj$mse.app, obj$mse.val))
      plot(x, obj$mse.app, xlab = xlab, ylab = blab, ylim = ylim, 
           type = "n")
      lines(x, obj$mse.app, lty = 3)
      lines(x, obj$mse.val, lty = 1)
      title(sub = paste(obj$xval, "-fold cross-validation", sep = ""),
            adj = 0)
      if(is.function(legendloc))
        legend(legendloc(1), c("Apparent", "Cross-validated"),
               lty = c(3, 1), bty = "n") else
      {
        par(usr=c(0,1,0,1))
        legend(legendloc[1],legendloc[2],
               c("Apparent", "Cross-validated"),
               lty = c(3, 1), bty = "n")
      }
    }
  if("dxy" %in% what)
    {
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
      else
        {
          par(usr=c(0,1,0,1))
          legend(legendloc[1],legendloc[2],
                 c("Apparent", "Cross-validated"),
                 lty = c(3, 1), bty = "n")
        }
    }
  invisible()
}
