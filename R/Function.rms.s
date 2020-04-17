Function.rms <- function(object, intercept=NULL,
                         digits=max(8,.Options$digits),
                         posterior.summary=c('mean', 'median', 'mode'), ...)
{
  posterior.summary <- match.arg(posterior.summary)
  
  oldopt <- options('digits')
  options(digits=digits)
  on.exit(options(oldopt))

  at   <- object$Design
  name <- at$name
  ac   <- at$assume.code
  p    <- length(name)
  nrp  <- num.intercepts(object)
  name.main <- name[ac!=9]  #non-intercepts
  pm        <- length(name.main)
  adj.to    <- Getlim(at, allow.null=TRUE, need.all=TRUE)$limits['Adjust to',]

  draws <- object$draws    # uses coef.rmsb if Bayesian
  Coef <- if(length(draws)) coef(object, stat=posterior.summary)
          else
            object$coef

  chr <- function(y, digits) if(is.factor(y) || is.character(y)) 
    paste('"',as.character(y),'"',sep='') else formatSep(y, digits)

  adj.to <- unlist(lapply(adj.to,chr,digits=digits))
  z <- paste('function(',paste(name.main,'=',adj.to,collapse=','), ') {', sep='')


  ##f$term.labels does not include strat
  TL <- attr(terms(object),"term.labels")
  ##Get inner transformations
  ##from <- c("asis","pol","lsp","rcs","catg","scored","strat","matrx","I")
  ##from <- paste(from,"(\\(.*\\))",sep="")
  from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
            'strat(*)','matrx(*)','I(*)')
  to   <- rep('*',9)

  ##trans <- paste("h(",translate(TL[ac!=9], from, "\\1"),")",sep="")  
  trans <- paste("h(",sedit(TL[ac!=9], from, to),")",sep="")
  ##change wrapping function to h()
  h <- function(x,...) deparse(substitute(x))
  for(i in (1:pm)) trans[i] <- eval(parse(text=trans[i]))
  j <- trans != name.main
  if(any(j)) z <- paste(z, paste(name.main[j],'<-',trans[j],collapse=';'),
                        ';',sep='')

  interaction <- at$interactions
  if(length(interaction) == 0) interaction <- 0
  
  parms <- at$parms
  
  Two.Way <- function(prm,Nam,nam.coef,cof,coef,f,varnames,at,digits)
    {
      i1 <- prm[1,1]; i2 <- prm[2,1]
      num.nl <- any(prm[1,-1] != 0)+any(prm[2,-1] != 0)
      ##If single factor with nonlinear terms, get it as second factor
      ##Otherwise, put factor with most # terms as second factor
      rev <- FALSE
      if((num.nl==1 & any(prm[1,-1] != 0)) ||
         (length(Nam[[i1]]) > length(Nam[[i2]])))
        { i1 <- i2; i2 <- prm[1,1]; rev <- TRUE }
      N1 <- Nam[[i1]]; N2 <- Nam[[i2]]
      n1 <- nam.coef[[i1]]; n2 <- nam.coef[[i2]]
      v <- ""
      for(j1 in 1:length(N1))
        {
          nam1 <- nam.coef[[i1]][j1]
          lN2 <- length(N2)
          cnam <- if(rev) paste(nam.coef[[i2]],"*",nam1) else
          paste(nam1, "*", nam.coef[[i2]])
          mnam <- match(cnam, names(cof), nomatch=0)
          act <- mnam[mnam>0]
          lN2.act <- length(act)
          ##Check if restricted interaction between a rcs and another nonlinear
          ##var, i.e. >1 2nd term possible, only 1 (linear) there, and at first 
          ##nonlinear term of rcs
          if(lN2.act==1 & lN2>1 & at$assume.code[i1]==4 & j1==2)
            {
              v <- paste(v,"+",N2[1],"*(",sep="")
              cnam <- paste(nam.coef[[if(rev)i2 else i1]][1], "*",
                            nam.coef[[if(rev)i1 else i2]][-1])
              vv <- attr(rcspline.restate(at$parms[[at$name[i1]]],
                                          c(0, coef[cnam]), 
                                          x=varnames[i1], digits=digits),
                         'function.text')
              v <- paste(v, vv, ')', sep='')
              break
            }
          else
            if(lN2.act==1)
              {
                vv <- paste(cof[act],"*",N1[j1],"*", N2[mnam>0], sep="")
                v <- paste(v, vv, sep='')
              }
            else
              if(lN2.act>0)
                {
                  vv <- paste("+",N1[j1],"*(",sep="")
                  v <- paste(v, vv, sep='')
                  
                  if(at$assume.code[i2]==4 & !any(mnam==0))
                    {
                      ##rcspline, interaction not restricted
                      vv <- attr(rcspline.restate(at$parms[[at$name[i2]]],
                                                  coef[act], 
                                                  x=varnames[i2],
                                                  digits=digits),
                                 'function.text')
                      v <- paste(v, vv, ')', sep='')
                    }
                  else
                    {
                      for(j2 in 1:lN2)
                        {
                          l <- mnam[j2]
                          if(l>0)
                            {	    #not a restricted-out nonlinear term
                              if(j2==1 && substring(cof[l],1,1)=="+")
                                cof[l] <- substring(cof[l],2)
                              vv <- paste(cof[l],"*",N2[j2],sep="")
                              v <- paste(v, vv, sep='')
                            }
                        }
                      v <- paste(v, ")", sep='')
                    }
                }
        }
      v
    }
  
Three.Way <- function(prm,Nam,nam.coef,cof,coef,f,at,digits)
  {
    i1 <- prm[1,1]; i2 <- prm[2,1]; i3 <- prm[3,1]
    N1 <- Nam[[i1]]; N2 <- Nam[[i2]]; N3 <- Nam[[i3]]
    v <- ""; l <- 0
    for(j3 in 1:length(N3))
      {
        for(j2 in 1:length(N2))
          {
            for(j1 in 1:length(N1))
              {
                l <- l+1
                v <- paste(v,cof[l], "*", N1[j1], "*", N2[j2],
                           "*", N3[j3], sep="")
              }
          }
      }
    v
  }
  

  if(nrp==1 | length(intercept))
    {
      cof <- if(! length(intercept)) formatSep(Coef[1], digits) else 
        formatSep(intercept, digits)
      z <- paste(z, cof, sep='')
    }
  
  Nam <- list();  nam.coef <- list()
  assig <- object$assign

  for(i in (1:p)) {
    ass <- ac[i]
    nam <- name[i]
    prm <- at$parms[[nam]]
    if(any(ass==c(5,7,8))) prm <- chr(at$parms[[nam]],digits=digits)
    
    k <- assig[[TL[i]]]
    coef <- Coef[k]
    nam.coef[[i]] <- names(coef)
    cof <- formatSep(coef,digits)
    cof <- ifelse(coef<=0, cof, paste("+", cof, sep=""))
    
    switch(ass,
           {
             nam <- name[i]; Nam[[i]] <- nam
             q <- paste(cof, '*', nam, sep="")
           },
           
           {
             q <- ""; pow <- 1:prm
             nams <- ifelse(pow==1,nam,paste(nam,"^",pow,"",sep=""))
             Nam[[i]] <- nams
             for(j in pow) q <- paste(q, cof[j], "*", nams[j], sep="")
           },
           
           {  
             q <- paste(cof[1], "*", nam, sep="")
             nams <- nam
             kn <- formatSep(-prm,digits)
             for(j in 1:length(prm)) {
               zz <- paste("pmax(", nam, if(prm[j]<0) "+" else NULL, 
                           if(prm[j]!=0) kn[j] else NULL, 
                           ",0)", sep="")
               nams <- c(nams, zz)
               q <- paste(q, cof[j+1], "*", zz, sep="")
             }
             Nam[[i]] <- nams
           },
           
           {
             q <- attr(rcspline.restate(prm, coef, x=nam, digits=digits),
                       'function.text')
             if(coef[1]>=0) q <- paste('+',q,sep='')
             nn <- nam
             for(j in 1:(length(prm)-2)) {
               nam <- paste(nam, "'", sep=""); nn <- c(nn, nam)
             }
             Nam[[i]] <- nn       #Two.Way only needs first name
                                        #for 2nd-order ia with 1 d.f. (restr ia)
                                        #Three.Way needs original design matrix
           } ,
           {
             nn <- paste('(',nam,'==',prm[-1],')',sep='')
             Nam[[i]] <- nn
             q <- ''
             for(j in 1:(length(prm)-1)) {
               vv <- paste(cof[j], nn[j], sep="*")
               q <- paste(q, vv, sep="")
             }
           },
           
           q <- '',
           
           {
             q <- paste(cof[1], "*", nam, sep="")
             nams <- nam
             for(j in 3:length(prm)) {
               zz <- prm[j]
               vv <- paste(cof[j-1], "*(", nam, "==", zz, ")", sep="")
               nams <- c(nams, zz)
               q <- paste(q, vv, sep="")
             }
             Nam[[i]] <- nams
           },
           ##Strat factor doesn't exist as main effect, but keep variable
           ##names and their lengths if they will appear in interactions later
           {
             ## was if(!length(Nam[[i]]) && any...
             if(any(interaction==i)) {
               nam.coef[[i]] <- paste(name[i], "=", prm[-1], sep="")
               Nam[[i]] <- prm[-1]
             }
             q <- "" },
           
           {  
             if(prm[3,1] == 0) 
               q <- Two.Way(prm,Nam,nam.coef,cof,coef,object,
                            name, at, digits)
             else q <- Three.Way(prm,Nam,nam.coef,cof,coef,
                                 object,at, digits)
             
           }, 
           {
             nam <- names(coef)
             q <- ""
             nam <- paste("(", nam, ")", sep="")
             Nam[[i]] <- nam
             for(j in 1:length(prm)) {
               vv <- paste(cof[j], '*', nam[j], sep="")
               q <- paste(q, vv, sep="")
             }
           }) 
    z <- paste(z, q, sep='')
  }
  z <- paste(z, '}')
  eval(parse(text=z))
}

Function.cph <-
  function(object, intercept=-object$center, ...)
  Function.rms(object, intercept=intercept, ...)

sascode <- function(object, file="", append=FALSE)
{
  chr <- function(y) if(is.factor(y) || is.character(y))
    paste('"',as.character(y),'"',sep='') else as.character(y)
  
  n <- names(object)[names(object)!='']
  for(i in n) if(file=='') cat(i,'=',chr(object[[i]]),';\n')
  else
    cat(i,'=',chr(object[[i]]),';\n',file=file, append=append|i>1)

  tf <- tempfile()
  dput(object, file=tf)
  object <- scan(file=tf, what='', sep='\n', quiet=TRUE)
  object <- paste(paste(object[3:(length(object)-1)],collapse='\n'),';',sep='')


  ##com <- 'sed -e "s/pmax/max/g" -e "s/pmin/min/g" -e "s/==/=/g" 
  ##-e "s/<-/=/g" -e "s/\\^/\*\*/g"'
  ##w <- sys(com, w)
  object <- sedit(object, c('pmax','pmin','==','<-','^'),
                  c('max','min','=','=','**'),
                  wild.literal=TRUE)
  if(file=='') cat(object, sep='\n')
  else
    cat(object, sep="\n", file=file, append=TRUE)
  invisible()
}

perlcode <- function(object) {
  group_translate <- function(expr) {
    result <- vector("list", length(expr) - 1)
    for (i in 2:length(expr)) {
      result[[i-1]] <- convert(expr[[i]])
    }
    paste(result, collapse=";\n  ")
  }

  simple_translate <- function(expr) {
    paste(convert(expr[[2]]), as.character(expr[[1]]), convert(expr[[3]]))
  }

  exp_translate <- function(expr) {
    expr[[1]] <- "**"
    simple_translate(expr)
  }

  pmax_pmin_translate <- function(expr) {
    result <- vector("list", length(expr) - 1)
    for (i in 2:length(expr)) {
      result[[i-1]] <- convert(expr[[i]])
    }
    name <- substr(as.character(expr[[1]]), 2, 4)
    paste(name, "((", paste(result, collapse=", "), "))", sep="")
  }

  equal_translate <- function(expr) {
    perlop <- if (is.character(expr[[2]]) || is.character(expr[[3]])) "eq" else "=="
    lhs <- convert(expr[[2]])
    rhs <- convert(expr[[3]])
    sprintf("(%s %s %s) ? 1 : 0", lhs, perlop, rhs)
  }

  parenthesis_translate <- function(expr) {
    sprintf("(%s)", convert(expr[[2]]))
  }

  assign_translate <- function(expr) {
    expr[[1]] <- "="
    simple_translate(expr)
  }

  log_translate <- function(expr) {
    paste("log(", convert(expr[[2]]), ")", sep="")
  }

  R_to_perl <- list(
    "{" = group_translate,
    "-" = simple_translate,
    "+" = simple_translate,
    "*" = simple_translate,
    "/" = simple_translate,
    "^" = exp_translate,
    "==" = equal_translate,
    "(" = parenthesis_translate,
    "<-" = assign_translate,
    "pmax" = pmax_pmin_translate,
    "pmin" = pmax_pmin_translate,
    "log" = log_translate
  )

  variable_translate <- function(v) {
    sprintf("$%s", gsub("\\.", "_", v))
  }

  convert <- function(expr) {
    if (length(expr) == 1) {
      x <- as.character(expr)
      if (typeof(expr) == "symbol") {
        variable_translate(x)
      } else {
        if (is.character(expr)) {
          sprintf('"%s"', x)
        }
        else {
          x
        }
      }
    }
    else {
      op <- as.character(expr[[1]])
      if (typeof(expr[[1]]) == "symbol" && op %in% names(R_to_perl)) {
        f <- R_to_perl[[op]]
        f(expr)
      } else {
        stop("don't know how to convert operator: ", op)
      }
    }
  }

  f <- object

  if (typeof(f) != "closure") {
    stop("argument must be a function")
  }

  fargs <- formals(f)
  fbody <- body(f)
  function_name <- as.character(match.call()[[2]])
  if (length(function_name) > 1) {
    function_name <- "f"
  }

  result <- list(sprintf("use List::Util 'max', 'min';\nsub %s {", function_name))
  for (i in 1:length(names(fargs))) {
    v <- names(fargs)[[i]]
    result <- c(result, sprintf("my %s = $_[%d];", variable_translate(v), i-1))
  }
  result <- c(result, convert(fbody))
  paste(paste(result, collapse="\n  "), "}", sep="\n")
}
