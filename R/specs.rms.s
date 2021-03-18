#Print description of specifications.  Can come from individual variable
#created by dx, complete design created by Design(), or complete design
#carried forward in fit

specs <- function(fit, ...) UseMethod('specs')

specs.rms <- function(fit, long=FALSE, ...)
{
  Call <- if(length(fit$call)) fit$call else
  if(length(attr(fit,'call'))) attr(fit,'call') else attr(fit, 'formula')

  tl <- attr(fit$terms, "term.labels")
  if(!length(tl)) tl <- attr(terms(formula(fit)), 'term.labels')
  ass    <- fit$assign
  strata <- levels(fit$strata)
  if(is.null(fit$assume)) {
    d <- fit$Design
    fit <- d
  }
  assume   <- fit$assume
  if(is.null(assume)) stop("fit does not have design information")
  parms    <- fit$parms
  name     <- fit$name
  lim      <- fit$limits
  ia.order <- fit$ia.order
  label    <- fit$label
  units    <- fit$units
  
  if(length(ass)) {
    if(names(ass)[1] %in% c("(Intercept)", "Intercept"))
      ass[[1]] <- NULL
    names(ass) <- name[assume != "strata"]
  }
  f <- length(assume)
  d <- matrix("", nrow=f, ncol=3)
  d[,1] <- assume
  iint  <- 0
  jfact <- 0
  trans <- rep("", f)
# Pick off inner transformation of variable. To complete, need to
# evaluate h function
# from <- c("asis","pol","lsp","rcs","catg","scored","strat","matrx","I")
# from <- paste(from,"(\\(.*\\))",sep="")
# tl <- translate(tl, from, "\\1")
# tl <- paste("h(",tl,")",sep="")

  from <- c('asis(*)', 'pol(*)', 'lsp(*)', 'rcs(*)', 'catg(*)', 'scored(*)',
            'strat(*)', 'matrx(*)', 'gTrans(*)', 'I(*)')
  to   <- rep('*', 10)

  tl <- paste("h(", sedit(tl, from, to), ")", sep="")
  ##change wrapping function to h()

  h <- function(x, ...) deparse(substitute(x))
  for(i in 1 : f)
    {
      if(assume[i] == "interaction") iint <- iint+1
      else
        {
          tr <- eval(parse(text = tl[i]))
          if(tr != name[i]) trans[i] <- tr
        }
      len <- if(assume[i] == "strata") 0 else length(ass[[name[i]]])
      d[i,3] <- as.character(len)
      parmi  <- parms[[name[i]]]
      if(d[i,1] == "transform") d[i,2] <- "function"
      else
        {
          if(length(parmi))
            {
              if(d[i,1] == "interaction")
                {
                  i1 <- parmi[1, -1] != 0
                  i2 <- parmi[2, -1] != 0
                  i3 <- parmi[3, -1] != 0
                  if(parmi[3,1] == 0)
                    {   #2nd order interaction
                      iao <- 1 * (any(i1) & !any(i2)) +
                        2 * (! any(i1) & any(i2)) +
                          3 * (any(i1) & any(i2) & ! any(i1 & i2)) +
                            4 * any(i1 & i2)
                      d[i,2] <- c("linear x linear - AB",
                                "nonlinear x linear - f(A)B",
                                "linear x nonlinear - Ag(B)",
                                "Af(B) + Bg(A)",
                                "f(A,B) - all cross-products")[iao+1]
                    }
                  else			#3rd order
                    d[i,2] <- paste(if(any(i1)) "nonlinear" else "linear", "x",
                                    if(any(i2)) "nonlinear" else "linear", "x",
                                    if(any(i3)) "nonlinear" else "linear")
                  if(ncol(parmi) == 1)  d[i,2] <- " "
                }
	
              else
                {
                  lab <- ""
                  if(assume[i] == 'gTrans')
                    parmi <- ''
                  for(z in parmi)
                    if(is.character(z)) lab <- paste(lab, z)
                    else
                      lab <- paste(lab, signif(z, 5))
                  d[i,2] <- lab
                }
            }
        }
    }
  collab <- c("Assumption", "Parameters", "d.f.")
  if(any(trans != ""))
    {
      collab <- c("Transformation", collab)
      d <- cbind(trans, d)
    }
  
  if(any(name != label))
    {
      collab <- c("Label", collab)
      d <- cbind(label, d)
    }
  if(length(units) && any(units != ''))
    {
      collab <- c('Units', collab)
      unitsb <- rep('', length(assume))
      unitsb[assume != 'interaction'] <- units
      d <- cbind(unitsb, d)
    }
  dimnames(d) <- list(name, collab)
  
  structure(list(call=Call, how.modeled=d, limits=if(long) lim,
                 strata=strata),
            class='specs.rms')
}

print.specs.rms <- function(x, ...) {
  dput(x$call)
  cat('\n')
  print(x$how.modeled, quote=FALSE)
  if(length(x$limits)) {cat('\n'); print(x$limits)}
  if(length(x$strata)) {
    cat("\n        Strata\n\n")
    print(x$strata,quote=FALSE)
  }
  invisible()
}

