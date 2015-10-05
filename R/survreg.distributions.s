# SCCS @(#)survreg.distributions.s	4.3 11/19/92
#
# Create the survreg.distributions object
#
# Infinite mean in log logistic courtesy of Victor Moreno
# SERC, Institut Catala d'Oncologia  (V.Moreno@ico.scs.es)  9Feb98

# survival package defines basic quantile function ignoring link
# Actual quantile function called Quantile here, for SV4 or R

survreg.auxinfo <- list(
exponential = list(
    survival = function(times, lp, parms) exp(-times/exp(lp)),
    hazard = function(times, lp, parms) exp(-lp),
    quantile = function(p) log(-log(p)),
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) -logb(1-q)*exp(lp)
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) exp(lp),
    
    latex = function(...) '\\exp(-t/\\exp(X\\beta))'
  ),
  
extreme = list(
    survival = function(times, lp, parms) { 
		exp(-exp((times-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms[1])   #14Jun97
		exp((times-lp)/scale)/scale
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(-logb(1-q))
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		lp-.57722*exp(parms)
		},
    latex = function(scale) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("\\exp[-\\exp(",z,")]")
		z
		}
    ),

weibull = list(
    survival = function(times, lp, parms) { 
		t.trans <- logb(times)
		names(t.trans) <- format(times)
		exp(-exp((t.trans-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		names(t.trans) <- format(times)
		scale <- exp(parms[1])   #14Jun97
		ifelse(times==0,exp(-lp/scale)/scale,
                        exp((t.trans-lp)/scale)*t.deriv/scale)
		},
    quantile = function(p) log(-log(p)),
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(-logb(1-q))
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms, transform) {
		names(parms) <- NULL
		exp(lp)*gamma(exp(parms)+1)
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("\\exp[-\\exp(",z,")]")
		z
		}
    ),
                    
logistic = list(
    survival = function(times, lp, parms) { 
		1/(1+exp((times-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms)
		1/scale/(1+exp(-(times-lp)/scale))
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(q/(1-q))
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale){
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("[1+\\exp(",z,")]^{-1}")
		z
		}
    ),

loglogistic = list(
    survival = function(times, lp, parms) { 
		1/(1+exp((logb(times)-lp)/exp(parms)))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		scale <- exp(parms)
		names(t.trans) <- format(times)
		t.deriv/scale/(1+exp(-(t.trans-lp)/scale))
		},
    quantile = qlogis,
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*logb(q/(1-q))
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		if(exp(parms)>1) rep(Inf,length(lp)) else
			   exp(lp)*pi*exp(parms)/sin(pi*exp(parms))
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("[1+\\exp(",z,")]^{-1}")
		z
		}),
    
gaussian = list(
    survival = function(times, lp, parms) pnorm(- (times-lp)/exp(parms)),
    hazard = function(times, lp, parms) {
		scale <- exp(parms)
		z <- (times-lp)/scale
		dnorm(z) / scale / pnorm(- z)
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*qnorm(q)
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z
		}
    ),

lognormal = list(
    survival = function(times, lp, parms) { 
		t.trans <- logb(times)
		names(t.trans) <- format(times)
		pnorm(- (t.trans-lp)/exp(parms))
		},
    hazard = function(times, lp, parms) {
		t.trans <- logb(times)
		t.deriv <- 1/times
		scale <- exp(parms)
		names(t.trans) <- format(times)
		z <- (t.trans-lp)/scale
		t.deriv * dnorm(z) / scale / pnorm(- z)
		},
    quantile = qnorm,
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms)*qnorm(q)
		names(q) <- format(q)
		drop(exp(outer(lp, q, FUN=f, parms=parms)))
		},
    mean = function(lp, parms) {
		names(parms) <- NULL
		exp(lp+exp(2*parms)/2)
		},
    latex = function(scale) {
		yvar <- "\\log(t)"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-\\Phi(",z,")")
		z
		}
    ),
  
t = list(
    survival = function(times, lp, parms) {
		scale <- exp(parms[1])
		df <- parms[2]
		pt(- (times-lp)/scale,df)
		},
    hazard = function(times, lp, parms) {
		scale <- exp(parms[1])
		df <- parms[2]
		z <- (times-lp)/scale
		dt(z,df) / scale / pt(- z,df)
		},
    Quantile = function(q=.5, lp, parms) {
		names(parms) <- NULL
		f <- function(lp, q, parms) lp + exp(parms[1])*qt(q, parms[2])
		names(q) <- format(q)
		drop(outer(lp, q, FUN=f, parms=parms))
		},
    mean = function(lp, parms) lp,
    latex = function(scale,df) {
		yvar <- "t"
		z <- if(scale==1) paste(yvar,"-X\\beta") else paste(
			"\\frac{", yvar, "-X\\beta}{",format(scale),"}",sep="")
		z <- paste("1-T_{",df,"}(",z,")", sep="")
		z
      }
  )
 )
