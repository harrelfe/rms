latexrms <- function(object,
		file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
		append=FALSE, which=1:p, varnames, columns=65, prefix=NULL, 
		inline=FALSE, before=if(inline)"" else "& &", intercept, 
		pretrans=TRUE, digits=.Options$digits)
{
  f    <- object	
  at   <- f$Design
  name <- at$name
  ac   <- at$assume.code
  p    <- length(name)
  nrp  <- num.intercepts(f)
  
  ## f$term.labels does not include strat
  TL   <- attr(terms(f),"term.labels")
  tl   <- TL
  
  ##Get inner transformations

  from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
            'strat(*)','matrx(*)','I(*)')
  to   <- rep('*',9)

  TLi <- paste("h(",sedit(TL, from, to),")",sep="")
  
  ## change wrapping function to h()

  h <- function(x,...) deparse(substitute(x))
  for(i in (1:p)[ac!=9]) TLi[i] <- eval(parse(text=TLi[i]))
  TLi <- ifelse(TLi==name | ac==1 | ac==9, "", TLi)
  anytr <- any(TLi!="")
  if(!missing(varnames))
    {
      if(length(varnames)!=sum(ac!=9)) stop("varnames is wrong length")
      vn <- name
      vn[ac!=9] <- varnames
      varnames <- vn
      tl <- sedit(tl, name, varnames)
      if(anytr) TLi <- sedit(TLi, name, varnames)
    }
  else
    varnames <- name
  lnam <- nchar(varnames)

  ## digits at end of name -> subscript, change font

  vnames <- sedit(varnames, '*$', '_{\\mit *}', test=all.digits)

  if(is.character(which))
    {
      wh <- charmatch(which, name, 0)
      if(any(wh==0))stop(paste("variable name not in model:",
               paste(which[wh==0], collapse=" ")))
    }

  interaction <- at$interactions
  if(length(interaction)==0) interaction <- 0
  
  parms <- at$parms

  ##If any interactions to be printed, make sure all main effects are included

  ia <- ac[which]==9
  if(length(which) < p & any(ia))
    {
      for(i in which[ia]) which <- c(which,parms[[name[i]]][,1])
      which <- which[which>0]
      which <- sort(unique(which))
    } 

  
  from <- c('sqrt(*)',  'log(',  'I(*)', '1/(*)',   'pmin(', 'pmax(')
  to   <- c('\\sqrt{*}','\\log(','[*]',  '(*)^{-1}','\\min(','\\max(')
  tl  <- sedit(tl, from, to)
  tl <- sedit(tl, varnames, vnames)
  ltl <- nchar(tl)
  tl <- paste("{\\rm ", tl, "}", sep="")
  if(anytr)
    {
      TLi <- sedit(TLi, from, to)
      TLi <- sedit(TLi, varnames, vnames)
      TLi <- ifelse(TLi=="", "", paste("{\\rm ", TLi, "}", sep=""))
    }
  
  varnames <- paste("{\\rm ", vnames, "}", sep="")
  
  Two.Way <- function(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,varnames,
                      lnam, at, digits=digits)
    {
      i1 <- prm[1,1]
      i2 <- prm[2,1]
      num.nl <- any(prm[1,-1] != 0)+any(prm[2,-1] != 0)
      ##If single factor with nonlinear terms, get it as second factor
      ##Otherwise, put factor with most # terms as second factor
      rev <- FALSE
      if((num.nl==1 & any(prm[1,-1] != 0)) ||
         (length(Nam[[i1]]) > length(Nam[[i2]])))
        {
          i1 <- i2
          i2 <- prm[1,1]
          rev <- TRUE
        }
      N1 <- Nam[[i1]]; N2 <- Nam[[i2]]
      n1 <- nam.coef[[i1]]; n2 <- nam.coef[[i2]]
      q <- NULL; cur <- ""; m <- 0
      for(j1 in 1:length(N1))
        {
          nam1 <- nam.coef[[i1]][j1]
          l1 <- lNam[[i1]][j1]
          lN2 <- length(N2)
          cnam <- if(rev) paste(nam.coef[[i2]], "*", nam1) else
          paste(nam1, "*", nam.coef[[i2]])
          mnam <- match(cnam, names(cof), nomatch=0)
          act <- mnam[mnam>0]
          lN2.act <- length(act)
          ##Check if restricted interaction between a rcs and another nonlinear
          ##var, i.e. >1 2nd term possible, only 1 (linear) there, and at first 
          ##nonlinear term of rcs
          
          if(lN2.act==1 & lN2>1 & at$assume.code[i1]==4 & j1==2)
            {
              if(cur!="")
                {
                  q <- c(q, cur)
                  m <- 0
                  cur <- ""
                }
              v <- paste("+",N2[1],"[",sep="")
              n <- lNam[[i2]][1]
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste(cur, v, sep="")
              m <- m+n
              cnam <- paste(nam.coef[[if(rev)i2 else i1]][1], "*",
                            nam.coef[[if(rev)i1 else i2]][-1])
              v <- rcspline.restate(at$parms[[at$name[i1]]], c(0, coef[cnam]), 
                                    x=varnames[i1],
                                    lx=lnam[i1], columns=columns, before="",
                                    after="",
                                    begin=cur, nbegin=m, digits=digits)
              m <- attr(v, "columns.used")+1   #+1 for "]"
              v <- attr(v, "latex")
              j <- length(v)
              if(j>1) q <- c(q, v[-j])
              cur <- paste(v[j], "]")
              break
            }
          else if(lN2.act==1)
            {
              v <- paste(cof[act],"\\:",N1[j1],"\\:\\times\\:",
                         N2[mnam>0], sep="")
              n <- l1+lNam[[i2]][mnam>0]+2
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste(cur, v, sep="")
              m <- m + n
            }
          else if(lN2.act>0)
            {
              if(cur!="")
                {
                  q <- c(q, cur)
                  m <- 0
                  cur <- ""
                }
              v <- paste("+",N1[j1],"[",sep="")
              n <- l1 + 1
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste(cur, v, sep="")
              m <- m + n
              
              if(at$assume.code[i2]==4 & !any(mnam==0))
                {
                  ##rcspline, interaction not restricted
                  v <- rcspline.restate(at$parms[[at$name[i2]]],
                                        coef[act], x=varnames[i2],
                                        lx=lnam[i2],
                                        columns=columns, before="",
                                        after="", 
                                        begin=cur, nbegin=m,
                                        digits=digits)
                  m <- attr(v, "columns.used") + 1   #1 for "]"
                  v <- attr(v, "latex")
                  j <- length(v)
                  if(j>1) q <- c(q, v[-j])
                  cur <- paste(v[j],"]")
                }
              
              else
                {
                  for(j2 in 1:lN2)
                    {
                      l <- mnam[j2]
                      if(l>0)
                        {	#not a restricted-out nonlinear term
                          if(j2==1 && substring(cof[l],1,1)=="+")
                            cof[l] <- substring(cof[l],2)
                          v <- paste(cof[l],"\\:",N2[j2],sep="")
                          n <- lcof[l]+lNam[[i2]][j2]
                          if(m + n > columns)
                            {
                              q <- c(q, cur)
                              cur <- ""
                              m <- 0
                            }
                          cur <- paste(cur, v, sep="")
                          m <- m + n
                        }
                    }
                  cur <- paste(cur, "]")
                }
            }
        }
      if(cur!="") q <- c(q, cur)
      attr(q, "columns.used") <- m
      q
    }
  
  Three.Way <- function(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,at)
    {
      i1 <- prm[1,1]; i2 <- prm[2,1]; i3 <- prm[3,1]
      N1 <- Nam[[i1]]; N2 <- Nam[[i2]]; N3 <- Nam[[i3]]
      q <- NULL
      cur <- ""
      m <- 0
      l <- 0
      for(j3 in 1:length(N3))
        {
          for(j2 in 1:length(N2))
            {
              for(j1 in 1:length(N1))
                {
                  l <- l + 1
                  v <- paste(cof[l], "\\:", N1[j1], "\\:\\times\\:", N2[j2],
                             "\\:\\times\\:", N3[j3], sep="")
                  n <- lcof[l] + lNam[[i1]][j1]+lNam[[i2]][j2] +
                    lNam[[i3]][j3] + 3
                  if(m + n > columns)
                    {
                      q <- c(q, cur)
                      cur <- ""
                      m <- 0
                    }
                  cur <- paste(cur, v, sep="")
                  m <- m + n
                }
            }
        }
      q <- c(q, cur)
      attr(q, "columns.used") <- m
      q
    }
  
  if(!inline)
    {
      tex <- "\\begin{eqnarray*}"
      if(length(prefix))
        tex <- c(tex,
                 paste("\\lefteqn{",prefix,"=}\\\\",sep=""))
    } else tex <- NULL
  
  cur <- ""
  cols <- 0
  Coef <- f$coef
  if((length(which)==p)&& (nrp==1 | !missing(intercept)))
    {
      cof <- if(missing(intercept)) format(Coef[1]) else format(intercept)
      cur <- cof
      cols <- nchar(cof)
    }
  
  anybrace <- anyplus <- FALSE
  Nam <- lNam <- nam.coef <- list()
  
  for(i in (1:p)[which])
    {
      ass <- ac[i]
      nam <- varnames[i]
      prm <- at$parms[[at$name[i]]]
      if(ass %in% c(5,7,8))
        {
          if(ass==7) prm <- format(prm)
          oprm <- prm
          lprm <- nchar(prm)
          z <- substring(prm,1,1)=="["
          u <- !z & ass==7
          prm <- sedit(prm, c(' ','&','%'), c('\\ ','\\&','\\%'))
          prm <- ifelse(z | u, prm, paste("{\\rm ", prm, "}", sep=""))
          prm <- ifelse(z,paste(nam,"\\in ",prm),prm)
          prm <- ifelse(u,paste(nam,"=",prm),prm)
          lprm <- lprm + (z | u)*(lnam[i]+1)
          prm <- paste("\\{", prm, "\\}", sep="")
          anybrace <- TRUE
        }
      if(ass != 8)
        {
          k <- f$assign[[TL[i]]]
          coef <- Coef[k]
          nam.coef[[i]] <- names(coef)
          cof <- format(coef)
          lcof <- nchar(cof)
          cof <- latexSN(cof)
          cof <- ifelse(coef<=0, cof, paste("+", cof, sep=""))
          cof.sp <- cof
          if(ass==2 | ass==10)
            {
              r <- grep("times",cof)
              r <- if(length(r)==0) 1:length(cof) else -r
              cof.sp[r] <- paste(cof.sp[r],"\\:",sep="")    
            }
          else
            if(length(grep("time",cof[1]))==0)
              cof.sp[1] <- paste(cof[1],"\\:",sep="")
          ## medium space between constant and variable names if constant
          ## does not end in 10^x
        }
      newline <- FALSE
      switch(ass,
             { # 1 - asis (linear)
               nam <- tl[i]
               Nam[[i]] <- nam
               lNam[[i]] <- ltl[i]
               q <- paste(cof.sp, nam, sep="")
               m <- ltl[i]+lcof
             },
             
             { # 2 - pol
               q <- ""
               m <- 0
               pow <- 1:prm
               nams <- ifelse(pow==1,nam,paste(nam,"^{",pow,"}",sep=""))
               Nam[[i]] <- nams; lNam[[i]] <- rep(lnam[i],prm)
               for(j in pow) q <- paste(q,cof.sp[j], nams[j], sep="")
               m <- prm*lnam[i]+sum(lcof)
             },
             
             { # 3 - lsp
               if(cols>0)
                 {
                   tex <- c(tex, cur)
                   cur <-""
                   cols <- 0
                 }
               anyplus <- TRUE
               q <- paste(cof.sp[1], nam, sep="")
               m <- lcof[1]+lnam[i]
               nams <- nam; lnams <- lnam[i]
               kn <- format(-prm)
               lkn <- nchar(kn)
               for(j in 1:length(prm))
                 {
                   z <- paste("(", nam, if(prm[j]<0) "+" else NULL, 
                              if(prm[j]!=0) kn[j] else NULL, ")_{+}",
                              sep="")
                 nams <- c(nams, z)
                   u <- lnam[i]+lkn[j]+2
                   lnams <- c(lnams,u)
                   v <- paste(cof[j+1], z, sep="")
                   n <- lcof[j+1]+u
                   if(m + n > columns)
                     {
                       cur <- paste(cur, q)
                       tex <- c(tex, cur)
                       cur <- ""
                       cols <- 0
                       q <- ""
                       m <- 0
                     }
                   q <- paste(q, v, sep="")
                   m <- m + n
                 }
               Nam[[i]] <- nams; lNam[[i]] <- lnams
             },
             
             
             { # 4 - rcs
               q <- rcspline.restate(prm, coef, x=nam, lx=lnam[i],
                                     columns=columns,
                                     before="",after="",digits=digits)
               anyplus <- TRUE
               m <- attr(q, "columns.used")
               nn <- nam; ln <- lnam[i]
               for(j in 1:(length(prm)-2))
                 {
                   nam <- paste(nam, "'", sep="")
                   nn <- c(nn, nam)
                   ln <- c(ln, lnam[i]+j)
                 }
               Nam[[i]] <- nn       #Two.Way only needs first name
               lNam[[i]] <- ln      #for 2nd-order ia with 1 d.f. (restr ia)
               ##Three.Way needs original design matrix
               q <- attr(q, "latex")
               if(substring(sedit(q[1]," ",""),1,1)!="-")
                 q[1] <- paste("+", q[1], sep="")
               j <- length(q)
               if(cur!="")
                 {
                   tex <- c(tex,cur)
                   cur <- ""
                   cols <- 0
                 }
               if(j>1)
                 {
                   tex <- c(tex, q[-j])
                   q <- q[j]
                 }
             } ,
             { # 5 - catg
               Nam[[i]] <- prm[-1]
               lNam[[i]] <- lprm[-1]
               if(cols>0)
                 {
                   tex <- c(tex,cur)
                   cur <- ""
                   cols <- 0
                 }
               q <- ""
               m <- 0
               for(j in 2:length(prm))
                 {
                   v <- paste(cof[j-1], prm[j], sep="")
                   n <- lcof[j-1]+lprm[j]
                   if(m + n > columns)
                     {
                       cur <- paste(cur,q)
                       tex <- c(tex, cur)
                       cur <- ""
                       cols <- 0
                       q <- ""
                       m <- 0
                     }
                   q <- paste(q, v, sep="")
                   m <- m + n
                 }
             },
             q <- "",
             
             { # 7 - scored
               if(cols>0)
                 {
                   tex <- c(tex,cur)
                   cur <- ""
                   cols <- 0
                 }
               q <- paste(cof.sp[1], nam, sep="")
               m <- nchar(q)
               nams <- nam
               lnams <- lnam[i]
               for(j in 3:length(prm))
                 {
                   z <- prm[j]
                   v <- paste(cof[j-1], z, sep="")
                   u <- lprm[j]+lnam[i]+3
                   n <- lcof[j-1]+u
                   nams <- c(nams, z)
                   lnams <- c(lnams,u)
                   if(m + n > columns)
                     {
                       cur <- paste(cur, q)
                       tex <- c(tex, cur)
                       cur <- ""
                       cols <- 0
                       q <- ""
                       m <- 0
                     }
                   q <- paste(q, v, sep="")
                   m <- m + n
                 }
               Nam[[i]] <- nams; lNam[[i]] <- lnams
             },
             ##Strat factor doesn't exist as main effect, but keep variable
             ##names and their lengths if they will appear in interactions later
             { # 8 - strat
               ## if(length(Nam[[i]])==0 && any(interaction==i)) 22Nov10
               if(any(interaction == i))
                 {
                   nam.coef[[i]] <- paste(name[i], "=", oprm[-1], sep="")
                   Nam[[i]] <- prm[-1]
                   lNam[[i]] <- lprm[-1]
                 }
               q <- ""
             },
             
             {
               if(prm[3,1]==0) 
                 q <- Two.Way(prm,Nam,nam.coef,lNam,cof,coef,f,columns,lcof,
                              varnames,lnam,at,digits=digits)
               else q <- Three.Way(prm,Nam,nam.coef,lNam,cof,coef,f,
                                   columns,lcof,at)
               m <- attr(q, "columns.used")
               j <- length(q)
               if(cur!="")
                 {
                   tex <- c(tex,cur)
                   cur <- ""
                   cols <- 0
                 }
               if(j>1)
                 {
                   tex <- c(tex,q[-j])
                   q <- q[j]
                 }
             }, 
             { # 10 - matrx
               nam <- names(coef)
               if(cols>0)
                 {
                   tex <- c(tex,cur)
                   cur <- ""
                   cols <- 0
                 }
               q <- ""
               m <- 0
               lnam <- nchar(nam)
               nam <- paste("{\\rm ", nam, "}", sep="")
               Nam[[i]] <- nam; lNam[[i]] <- lnam
               for(j in 1:length(prm))
                 {
                   v <- paste(cof.sp[j], nam[j], sep="")
                 n <- lcof[j]+lnam[j]
                   if(m + n > columns)
                     {
                       cur <- paste(cur, q)
                       tex <- c(tex, cur)
                       cur <- ""
                       cols <- 0
                       q <- ""
                       m <- 0
                     }
                   q <- paste(q, v, sep="")
                   m <- m + n
                 }
             }
             ) 
    
     
      if(length(q) && q!="")
        {
          if(cols+m > columns)
            {
              tex <- c(tex, cur)
              cur <- ""
              cols <- 0
            }
          
          cur <- paste(cur, q)
          cols <- cols + m
        }
    }
  
  if(cur!="") tex <- c(tex, cur)
  
  if(inline)
    {
      cat(tex, sep="\n", file=file, append=append)
      return(structure(list(file=file,style=NULL), class='latex'))
    }
  
  tex <- c(tex,"\\end{eqnarray*}")
  tex <- ifelse(substring(tex,1,1)=="\\",tex,paste(before,tex,"\\\\"))
  
  if(anybrace | anyplus)
    {
      s <- if(length(which)==p) "and $" else "where $"
      if(anybrace)
        s <- paste(s,"\\{c\\}=1 {\\rm\\ if\\ subject\\ is\\ in\\ group\\ } c, \\ 0 {\\rm\\ otherwise}")
      if(anybrace & anyplus) s <- paste(s, ";\\ ")
      if(anyplus)
        s <- paste(s, "(x)_{+}=x {\\rm\\ if\\ } x>0, \\ 0 {\\rm\\ otherwise}")
      s <- paste(s, "$.")
      tex <- c(tex, s)
    }
  
  if(anytr & pretrans)
    {
      i <- TLi!=""
      if(sum(i)==1) tr <- paste("$",varnames[i],
              "$ is pre--transformed as $",TLi[i],"$.",sep="")
      else
        {
          tr <- c("\\vspace{0.5ex}\\begin{center}{\\bf Pre--Transformations}\\\\",
                  "\\vspace{1.5ex}\\begin{tabular}{|l|l|} \\hline",
                  "\\multicolumn{1}{|c|}{Variable} & \\multicolumn{1}{c|}{Transformation} \\\\ \\hline",
                  paste("$",varnames[i],"$ & $",TLi[i],"$ \\\\",sep=""),
                  "\\hline", "\\end{tabular}\\end{center}")
        }
      tex <- c(tex, tr)
    }
  
  cat(tex, sep="\n", file=file, append=append)
  structure(list(file=file, style=NULL),class='latex')
}
