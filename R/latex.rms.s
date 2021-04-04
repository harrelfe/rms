latexrms <-
  function(object,
           file="",
           append=FALSE, which=1 : p, varnames, columns=65, prefix=NULL, 
           inline=FALSE, before=if(inline) "" else "& &", after="",
           intercept, pretrans=TRUE, digits=.Options$digits, size='')
{
  html <- prType() == 'html'
  ## Break character for non-math mode:
  brchar <- if(html) '<br>' else '\\\\'
  
  f    <- object	
  at   <- f$Design
  name <- at$name
  ac   <- at$assume.code
  Tex  <- at$tex
  p    <- length(name)
  nrp  <- num.intercepts(f)
  
  ## f$term.labels does not include strat
  TL   <- attr(terms(f), "term.labels")
  tl   <- TL
  
  ##Get inner transformations

  from <- c('asis(*)','pol(*)','lsp(*)','rcs(*)','catg(*)','scored(*)',
            'strat(*)','matrx(*)','gTrans(*)','I(*)')
  to   <- rep('*',10)

  TLi <- paste0("h(",sedit(TL, from, to),")")
  
  ## change wrapping function to h()

  h <- function(x,...) deparse(substitute(x))
  for(i in (1 : p)[ac != 9]) TLi[i] <- eval(parse(text=TLi[i]))
  TLi <- ifelse(TLi == name | ac == 1 | ac == 9, "", TLi)
  anytr <- any(TLi != "")
  if(! missing(varnames)) {
    if(length(varnames) != sum(ac != 9)) stop("varnames is wrong length")
    vn <- name
    vn[ac != 9] <- varnames
    varnames <- vn
    tl <- sedit(tl, name, varnames, wild.literal=TRUE)
    if(anytr) TLi <- sedit(TLi, name, varnames, wild.literal=TRUE)
  }
  else
    varnames <- name
  lnam <- nchar(varnames)

  ## digits at end of name -> subscript, change font
  ## used to be {\\mit *}

  vnames <- sedit(varnames, '*$', '_{*}', test=all.digits)

  if(is.character(which))
    {
      wh <- charmatch(which, name, 0)
      if(any(wh == 0))stop(paste("variable name not in model:",
               paste(which[wh == 0], collapse=" ")))
    }

  interaction <- at$interactions
  if(length(interaction) == 0) interaction <- 0
  
  parms <- at$parms

  ##If any interactions to be printed, make sure all main effects are included

  ia <- ac[which] == 9
  if(length(which) < p & any(ia))
    {
      for(i in which[ia]) which <- c(which,parms[[name[i]]][,1])
      which <- which[which>0]
      which <- sort(unique(which))
    } 

  
  from <- c('sqrt(*)',  'log(',  'I(*)', '1/(*)',   'pmin(', 'pmax(')
  to   <- c('\\sqrt{*}','\\log(','[*]',  '(*)^{-1}','\\min(','\\max(')
  tl  <- sedit(tl, from, to)
  tl <- sedit(tl, varnames, vnames, wild.literal=TRUE)
  ltl <- nchar(tl)
  tl <- paste0("\\mathrm{", tl, "}")
  if(anytr)
    {
      TLi <- sedit(TLi, from, to)
      TLi <- sedit(TLi, varnames, vnames, wild.literal=TRUE)
      TLi <- ifelse(TLi == "", "", paste0("\\mathrm{", TLi, "}"))
    }
  
  varnames <- paste0("\\mathrm{", vnames, "}")
  
  Two.Way <- function(prm, Nam, nam.coef, lNam, cof, coef, f,
                      columns, lcof, varnames,
                      lnam, at, digits=digits)
    {
      i1 <- prm[1, 1]
      i2 <- prm[2, 1]
      num.nl <- any(prm[1, -1] != 0) + any(prm[2, -1] != 0)
      ##If single factor with nonlinear terms, get it as second factor
      ##Otherwise, put factor with most # terms as second factor
      rev <- FALSE
      if((num.nl == 1 & any(prm[1, -1] != 0)) ||
         (length(Nam[[i1]]) > length(Nam[[i2]])))
        {
          i1  <- i2
          i2  <- prm[1,1]
          rev <- TRUE
        }
      N1 <- Nam[[i1]];      N2 <- Nam[[i2]]
      n1 <- nam.coef[[i1]]; n2 <- nam.coef[[i2]]
      q <- NULL; cur <- ""; m <- 0
      for(j1 in 1 : length(N1))
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
          
          if(lN2.act == 1 & lN2>1 & at$assume.code[i1] == 4 & j1 == 2)
            {
              if(cur != "")
                {
                  q <- c(q, cur)
                  m <- 0
                  cur <- ""
                }
              v <- paste0("+", N2[1], "[")
              n <- lNam[[i2]][1]
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste0(cur, v)
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
          else if(lN2.act == 1)
            {
              v <- paste0(cof[act],"\\:",N1[j1],"\\:\\times\\:",
                          N2[mnam>0])
              n <- l1+lNam[[i2]][mnam > 0] + 2
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste0(cur, v)
              m <- m + n
            }
          else if(lN2.act>0)
            {
              if(cur != "")
                {
                  q <- c(q, cur)
                  m <- 0
                  cur <- ""
                }
              v <- paste0("+", N1[j1], "[")
              n <- l1 + 1
              if(m + n > columns)
                {
                  q <- c(q, cur)
                  cur <- ""
                  m <- 0
                }
              cur <- paste0(cur, v)
              m <- m + n
              
              if(at$assume.code[i2] == 4 & ! any(mnam == 0))
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
                  for(j2 in 1 : lN2)
                    {
                      l <- mnam[j2]
                      if(l > 0)
                        {	#not a restricted-out nonlinear term
                          if(j2 == 1 && substring(cof[l],1,1) == "+")
                            cof[l] <- substring(cof[l],2)
                          v <- paste0(cof[l], "\\:", N2[j2])
                          n <- lcof[l] + lNam[[i2]][j2]
                          if(m + n > columns)
                            {
                              q <- c(q, cur)
                              cur <- ""
                              m <- 0
                            }
                          cur <- paste0(cur, v)
                          m <- m + n
                        }
                    }
                  cur <- paste(cur, "]")
                }
            }
        }
      if(cur != "") q <- c(q, cur)
      attr(q, "columns.used") <- m
      q
    }
  
  Three.Way <- function(prm, Nam, nam.coef, lNam, cof, coef, f,
                        columns, lcof, at)
    {
      i1 <- prm[1,1];  i2 <- prm[2,1];  i3 <- prm[3,1]
      N1 <- Nam[[i1]]; N2 <- Nam[[i2]]; N3 <- Nam[[i3]]
      q <- NULL
      cur <- ""
      m <- 0
      l <- 0
      for(j3 in 1 : length(N3))
        {
          for(j2 in 1 : length(N2))
            {
              for(j1 in 1 : length(N1))
                {
                  l <- l + 1
                  v <- paste0(cof[l], "\\:", N1[j1], "\\:\\times\\:", N2[j2],
                              "\\:\\times\\:", N3[j3])
                  n <- lcof[l] + lNam[[i1]][j1] + lNam[[i2]][j2] +
                    lNam[[i3]][j3] + 3
                  if(m + n > columns)
                    {
                      q <- c(q, cur)
                      cur <- ""
                      m <- 0
                    }
                  cur <- paste0(cur, v)
                  m <- m + n
                }
            }
        }
      q <- c(q, cur)
      attr(q, "columns.used") <- m
      q
    }
  
  if(! inline)
    {
      tex <- "\\begin{eqnarray*}"
      if(size != '') tex <- c(tex, paste0('\\', size))
      if(length(prefix))
        tex <- c(tex,
                 if(html) paste0(prefix, '= & & \\\\') else
                 paste0("\\lefteqn{", prefix, "=}\\\\"))
    } else tex <- NULL
  
  cur <- ""
  cols <- 0
  Coef <- f$coef
  if((length(which) == p)&& (nrp == 1 | ! missing(intercept)))
    {
      cof <- if(missing(intercept))
        format(Coef[1], digits=digits) else format(intercept, digits=digits)
      cur <- cof
      cols <- nchar(cof)
    }
  
  anyivar <- anyplus <- FALSE   # anyivar = any indicator variable
  Nam <- lNam <- nam.coef <- list()
  
  for(i in (1 : p)[which])
    {
      ass <- ac[i]
      nam <- varnames[i]
      prm <- at$parms[[at$name[i]]]
      if(ass %in% c(5, 7, 8))
        {
          if(ass == 7) prm <- format(prm)
          oprm <- prm
          lprm <- nchar(prm)
          z <- substring(prm, 1, 1) == "["
          u <- ! z & ass == 7
          prm <- sedit(prm, c(' ','&','%'), c('\\ ','\\&','\\%'))
          prm <- ifelse(z | u, prm, paste0("\\mathrm{", prm, "}"))
          prm <- ifelse(z, paste(nam, "\\in ", prm), prm)
          prm <- ifelse(u, paste(nam, "=", prm), prm)
          lprm <- lprm + (z | u) * (lnam[i] + 1)
          prm <- paste0("[", prm, "]")
          anyivar <- TRUE
        }
      if(ass != 8)
        {
          k <- f$assign[[TL[i]]]
          coef <- Coef[k]
          nam.coef[[i]] <- names(coef)
          cof <- formatSep(coef, digits=digits)
          lcof <- nchar(cof)
          cof <- latexSN(cof)
          cof <- ifelse(coef<=0, cof, paste0("+", cof))
          cof.sp <- cof
          if(ass == 2 | ass == 10)
            {
              r <- grep("times", cof)
              r <- if(length(r) == 0) 1 : length(cof) else -r
              cof.sp[r] <- paste0(cof.sp[r], "\\:")    
            }
          else
            if(length(grep("time",cof[1])) == 0)
              cof.sp[1] <- paste0(cof[1], "\\:")
          ## medium space between constant and variable names if constant
          ## does not end in 10^x
        }
      newline <- FALSE
      switch(ass,
             { # 1 - asis (linear)
               nam       <- tl[i]
               Nam[[i]]  <- nam
               lNam[[i]] <- ltl[i]
               q <- paste0(cof.sp, nam)
               m <- ltl[i] + lcof
             },
             
             { # 2 - pol
               q   <- ""
               m   <- 0
               pow <- 1 : prm
               nams <- ifelse(pow == 1,nam, paste0(nam, "^{", pow, "}"))
               Nam[[i]] <- nams; lNam[[i]] <- rep(lnam[i],prm)
               for(j in pow) q <- paste0(q,cof.sp[j], nams[j])
               m <- prm * lnam[i] + sum(lcof)
             },
             
             { # 3 - lsp
               if(cols > 0)
                 {
                   tex <- c(tex, cur)
                   cur <-""
                   cols <- 0
                 }
               anyplus <- TRUE
               q <- paste0(cof.sp[1], nam)
               m <- lcof[1] + lnam[i]
               nams <- nam; lnams <- lnam[i]
               kn <- format(-prm)
               lkn <- nchar(kn)
               for(j in 1 : length(prm))
                 {
                   z <- paste0("(", nam, if(prm[j] < 0) "+" else NULL, 
                              if(prm[j] != 0) kn[j] else NULL, ")_{+}")
                 nams <- c(nams, z)
                   u <- lnam[i] + lkn[j] + 2
                   lnams <- c(lnams, u)
                   v <- paste0(cof[j + 1], z)
                   n <- lcof[j + 1] + u
                   if(m + n > columns)
                     {
                       cur  <- paste(cur, q)
                       tex  <- c(tex, cur)
                       cur  <- ""
                       cols <- 0
                       q    <- ""
                       m    <- 0
                     }
                   q <- paste0(q, v)
                   m <- m + n
                 }
               Nam[[i]] <- nams; lNam[[i]] <- lnams
             },
             
             
             { # 4 - rcs
               q <- rcspline.restate(prm, coef, x=nam, lx=lnam[i],
                                     columns=columns,
                                     before="", after="", digits=digits)
               anyplus <- TRUE
               m <- attr(q, "columns.used")
               nn <- nam; ln <- lnam[i]
               for(j in 1 : (length(prm) - 2))
                 {
                   nam <- paste0(nam, "'")
                   nn <- c(nn, nam)
                   ln <- c(ln, lnam[i] + j)
                 }
               Nam[[i]]  <- nn      #Two.Way only needs first name
               lNam[[i]] <- ln      #for 2nd-order ia with 1 d.f. (restr ia)
               ##Three.Way needs original design matrix
               q <- attr(q, "latex")
               if(substring(sedit(q[1], " ", ""), 1, 1) != "-")
                 q[1] <- paste0("+", q[1])
               j <- length(q)
               if(cur != "")
                 {
                   tex  <- c(tex,cur)
                   cur  <- ""
                   cols <- 0
                 }
               if(j > 1)
                 {
                   tex <- c(tex, q[-j])
                   q   <- q[j]
                 }
             } ,
             { # 5 - catg
               Nam[[i]]  <- prm[-1]
               lNam[[i]] <- lprm[-1]
               if(cols > 0)
                 {
                   tex  <- c(tex,cur)
                   cur  <- ""
                   cols <- 0
                 }
               q <- ""
               m <- 0
               for(j in 2 : length(prm))
                 {
                   v <- paste0(cof[j - 1], prm[j])
                   n <- lcof[j - 1] + lprm[j]
                   if(m + n > columns)
                     {
                       cur  <- paste(cur, q)
                       tex  <- c(tex, cur)
                       cur  <- ""
                       cols <- 0
                       q    <- ""
                       m    <- 0
                     }
                   q <- paste0(q, v)
                   m <- m + n
                 }
             },
             q <- "",
             
             { # 7 - scored
               if(cols > 0)
                 {
                   tex <- c(tex, cur)
                   cur <- ""
                   cols <- 0
                 }
               q <- paste0(cof.sp[1], nam)
               m <- nchar(q)
               nams <- nam
               lnams <- lnam[i]
               for(j in 3 : length(prm))
                 {
                   z <- prm[j]
                   v <- paste0(cof[j - 1], z)
                   u <- lprm[j] + lnam[i] + 3
                   n <- lcof[j - 1] + u
                   nams  <- c(nams, z)
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
                   q <- paste0(q, v)
                   m <- m + n
                 }
               Nam[[i]] <- nams; lNam[[i]] <- lnams
             },
             ##Strat factor doesn't exist as main effect, but keep variable
             ##names and their lengths if they will appear in interactions later
             { # 8 - strat
               ## if(length(Nam[[i]]) == 0 && any(interaction == i)) 22Nov10
               if(any(interaction == i))
                 {
                   nam.coef[[i]] <- paste0(name[i], "=", oprm[-1])
                   Nam[[i]]  <- prm[-1]
                   lNam[[i]] <- lprm[-1]
                 }
               q <- ""
             },
             
             {
               if(prm[3,1] == 0) 
                 q <- Two.Way(prm, Nam, nam.coef, lNam, cof, coef, f,
                              columns, lcof,
                              varnames, lnam, at, digits=digits)
               else q <- Three.Way(prm, Nam, nam.coef, lNam, cof, coef, f,
                                   columns, lcof, at)
               m <- attr(q, "columns.used")
               j <- length(q)
               if(cur != "")
                 {
                   tex  <- c(tex,cur)
                   cur  <- ""
                   cols <- 0
                 }
               if(j > 1)
                 {
                   tex <- c(tex,q[-j])
                   q   <- q[j]
                 }
             }, 
             { # 10 - matrx
               nam <- names(coef)
               if(cols > 0)
                 {
                   tex  <- c(tex,cur)
                   cur  <- ""
                   cols <- 0
                 }
               q <- ""
               m <- 0
               lnam <- nchar(nam)
               nam <- paste0("\\mathrm{", nam, "}")
               Nam[[i]] <- nam; lNam[[i]] <- lnam
               for(j in 1 : length(prm))
                 {
                   v <- paste0(cof.sp[j], nam[j])
                 n <- lcof[j] + lnam[j]
                   if(m + n > columns)
                     {
                       cur  <- paste(cur, q)
                       tex  <- c(tex, cur)
                       cur  <- ""
                       cols <- 0
                       q    <- ""
                       m    <- 0
                     }
                   q <- paste0(q, v)
                   m <- m + n
                 }
             },

             { # 11 - gTrans
               if(! length(Tex))
                 stop(ta <-
                        paste('no tex attribute for gTrans variable', name[i]))
               tx   <- Tex[[name[i]]]
               if(! length(tx)) stop(z)
               tx   <- eval(parse(text=tx))
               nams <- tx(nam)
               q <- ""
               m <- 0
               lx <- length(nams)
               Nam[[i]] <- nams; lNam[[i]] <- rep(lnam[i], lx)
               for(j in 1 : lx) q <- paste0(q, cof.sp[j], nams[j])
               m <- lx * lnam[i] + sum(lcof)
             }
             ) 
    
     
      if(length(q) && q != "")
        {
          if(cols + m > columns)
            {
              tex <- c(tex, cur)
              cur <- ""
              cols <- 0
            }
          
          cur <- paste(cur, q)
          cols <- cols + m
        }
    }
  
  if(cur != "") tex <- c(tex, cur)

  if(inline) {
    if(before != '') tex <- c(before, tex)
    if(size != '')   tex <- c(paste0('{\\', size), tex)
    if(after  != '') tex <- c(tex, after)
    if(size != '')   tex <- c(tex, '}')
    if(html) return(htmltools::HTML(paste0(tex, '\n')))
    cat(tex, sep="\n", file=file, append=append)
    return(structure(list(file=file,style=NULL), class='latex'))
  }
  
  tex <- c(tex, "\\end{eqnarray*}")

  tex <- ifelse(tex == paste0(prefix, '= & & \\\\') |
                substring(tex,1,1) == "\\", tex,
                paste(before, tex, "\\\\"))
  
  if(anyivar | anyplus) {
    s <- if(length(which) == p) "and " else "where "
    if(anyivar)
      s <- paste0(s, "\\([c]=1\\) if subject is in group \\(c\\), 0 otherwise")
    ## Had trouble with Rmarkdown recognizing math mode with $...$
    if(anyivar && anyplus) s <- paste0(s, '; ')
    if(anyplus)
      s <- paste0(s, "\\((x)_{+}=x\\) if \\(x > 0\\), 0 otherwise", brchar)
    tex <- c(tex, s)
  }
  
  if(anytr & pretrans) {
    i <- TLi != ""
    if(sum(i) == 1) tr <- paste0("\\(", varnames[i],
                                 "\\) is pre--transformed as \\(",
                                 TLi[i], "\\).")
    else {
      tr <- if(html) {
              z <- cbind(Variable=paste0('\\(', varnames, '\\)'),
                         Transformation=paste0('\\(', TLi, '\\)'))
              as.character(htmlTable::htmlTable(z, caption='Pre-transformations',
                                                css.cell='min-width: 9em;',
                                                align='|l|l|',
                                                align.header='|c|c|',
                                                escape.html=FALSE))
#                           sep='\n')
            }
            else
              c("\\vspace{0.5ex}\\begin{center}{\\bf Pre--Transformations}\\\\",
                "\\vspace{1.5ex}\\begin{tabular}{|l|l|} \\hline",
                "\\multicolumn{1}{|c|}{Variable} & \\multicolumn{1}{c|}{Transformation} \\\\ \\hline",
                paste0("\\(",varnames[i],"\\) & \\(",TLi[i],"\\) \\\\"),
                "\\hline", "\\end{tabular}\\end{center}")
    }
    tex <- c(tex, tr)
  }
  if(html) return(htmltools::HTML(paste0(tex, '\n')))
  cat(tex, sep="\n", file=file, append=append)
  structure(list(file=file, style=NULL), class='latex')
}
