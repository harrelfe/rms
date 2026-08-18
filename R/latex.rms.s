## -----------------------------------------------------------------------
## latexrms, with changes concentrated in four places, each marked below.
## Two.Way, Three.Way, and the entire term-building switch (assume.code
## cases 1-11) are byte-for-byte identical to the uploaded source --
## nothing there needed to change, since none of it emits LaTeX syntax
## directly; it builds raw term/coefficient content that the
## before/after wrapping (at the end) and the array scaffolding (at the
## start) apply uniformly to.
##
## No prType()-specific branches were added anywhere in these changes --
## per the conclusion that plain LaTeX math, made texmath-safe, renders
## correctly under both prType='latex' (real LaTeX/PDF) and
## prType='typst' (via Pandoc's texmath conversion) from the SAME code
## path. The existing html-specific pre-transformations-table branch
## near the end (checking `html <- prType()=='html'`) is untouched --
## that's a separate, already-working feature, out of scope here.
##
## CHANGE 1 (signature): before=if(inline) "" else "& &" -> "&". The
## array is now 2 columns (see Change 2), so continuation rows need one
## leading empty column, not two.
##
## CHANGE 2 (array opening + prefix row): \begin{array} -> explicit
## \begin{array}{ll} -- confirmed by direct testing that texmath's array
## reader requires a column-spec argument and errors ("expecting ... {")
## on a bare \begin{array}. Left-aligned per explicit preference (a
## right-aligned version was tried, to line up the trailing = in
## "X\hat{\beta}=", "[c]=", and "(x)_+=", but reverted). \lefteqn{prefix=}\\
## -> prefix= & \\ -- \lefteqn{} is an eqnarray-package macro texmath's
## LaTeX reader doesn't recognize (confirmed by testing); the 2-column
## "prefix in column 1, empty column 2" row achieves the identical
## hanging-indent visual layout using only standard array semantics, no
## package-specific macro needed. The [c]=/(x)_+= definition rows (see
## CHANGE 5) use this same "label in column 1" structure -- an earlier
## version of CHANGE 5 incorrectly put them in column 2 instead, via the
## continuation-row 'before' prefix, which was the actual bug behind
## them appearing misaligned with the rest of the array.
##
## CHANGE 3 (inline branch): the \\\\-joining of tex elements now has a
## {} guard after every join, not just implicitly at the end -- same
## fix as changes 4/5 below, applied here too since it's the same
## silent-content-swallowing risk (a joined element starting with '['
## right after a bare \\ gets misparsed as \\[...], an optional
## LaTeX spacing argument, and silently disappears with no error).
##
## CHANGE 4 (closing section): tex == paste0(prefix, '= & & \\\\') check
## updated to match the new 2-column prefix row ('= & \\\\', not
## '= & & \\\\'). The substring(tex,1,1)=="\\" guard is kept as a
## general safety net even though \end{array} no longer passes through
## this ifelse() at all (see Change 5) -- harmless to keep, protects any
## future/edge-case row that happens to start with a literal backslash.
##
## CHANGE 5 (closing section): the anyivar/anyplus trailing definitions
## ([c]=1 if..., (x)_+=x if...) are now merged into the SAME array as
## additional rows, instead of being appended afterward as separate $$
## blocks -- confirmed empirically to produce much better, consistent
## spacing. The last coefficient/term row's normal single "\\\\"
## terminator gets extended to "\\\\{} \\\\{}" (an extra guarded blank
## row) to visually separate the two sections within the one array.
## \end{array} now closes the array AFTER these definition rows are
## appended, not before.
## -----------------------------------------------------------------------
latexrms <-
  function(object,
           file="",
           append=FALSE, which=1 : p, varnames, columns=65, prefix=NULL, 
           inline=FALSE, before=if(inline) "" else "&", after="",   ## CHANGE 1
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
      tex <- '\\begin{array}{ll}'   ## CHANGE 2: explicit 2-column spec,
                                    ## left-aligned
      if(size != '') tex <- c(tex, paste0('\\', size))
      if(length(prefix))
        tex <- c(tex,
                 paste0(prefix, "= & \\\\"))
                 ## CHANGE 2: was \lefteqn{prefix=}\\ -- see header note
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
          if(ass != 8) anyivar <- TRUE
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
    tex <- paste(tex, collapse='\\\\{}')   ## CHANGE 3: {} guard on every join
    tex <- c('\\begin{array}{l}', tex, '\\end{array}')
    if(before != '') tex <- c(before, '', tex)
    if(size != '')   tex <- c(paste0('{\\', size), tex)
    if(after  != '') tex <- c(tex, after)
    if(size != '')   tex <- c(tex, '}')
    if(file == '') return(tex)
    cat(tex, sep='\n', file=file, append=append)
    return(invisible())
  }
  
  tex <- ifelse(tex == paste0(prefix, '= & \\\\') |    ## CHANGE 4: '& \\\\' not '& & \\\\'
                substring(tex,1,1) == "\\", tex,
                paste(before, tex, "\\\\"))

  ## CHANGE 5: merge trailing definitions into the same array, with a
  ## guarded blank-row separator, instead of appending as separate $$
  ## blocks -- see header note.
  if(anyivar | anyplus) {
    tex[length(tex)] <- paste0(tex[length(tex)], "{} \\\\{}")

    ## Label goes in column 1 (matching the X\hat{\beta}= prefix row's
    ## own structure), NOT appended after 'before' into column 2 -- that
    ## was the actual bug: "before" is the empty-column-1 continuation
    ## prefix meant for ordinary coefficient rows, not appropriate here
    ## since these rows have their own label, the same as the prefix row.
    defs <- character(0)
    if(anyivar)
      defs <- c(defs,
                "[c]= & 1~\\mathrm{if~subject~is~in~group}~c,~0~\\mathrm{otherwise} \\\\")
    if(anyplus)
      defs <- c(defs,
                "(x)_{+}= & x~\\mathrm{if}~x > 0,~0~\\mathrm{otherwise} \\\\")
    tex <- c(tex, defs)
  }

  tex <- c(tex, '\\end{array}')  # was eqnarray*

  ## NEW: wrap the array in $$ ... $$ so it's recognized as genuine math
  ## -- by real LaTeX (which needs math mode for \hat{}/_{}/^{} to mean
  ## anything at all) and by Pandoc/texmath under Typst (which only
  ## translates content inside math delimiters). This was missing
  ## entirely before -- latexrms never added it, and none of the
  ## calling latex.X functions did either, so the array was always
  ## reaching the compiler as bare text outside math mode. Applied here
  ## by position (first/last element), after the before/after ifelse()
  ## step above, specifically to avoid interfering with that step's
  ## existing substring(tex,1,1)=="\\" protection -- prepending $$
  ## earlier would make the array-opening line start with '$' instead
  ## of '\\', failing that exclusion check and getting an unwanted
  ## before/after wrapping applied to it too. Only applied for the
  ## non-inline case -- inline mode is presumably meant to be embedded
  ## within a math context the caller already has open, though that
  ## assumption hasn't been specifically confirmed.
  if(! inline) {
    tex[1]            <- paste0('$$', tex[1])
    tex[length(tex)]  <- paste0(tex[length(tex)], '$$')
  }

  if(anytr & pretrans) {
    i <- TLi != ""
    if(sum(i) == 1) tr <- paste0("$", varnames[i],
                                 "$ is pre--transformed as $",
                                 TLi[i], "$.")
    else {
      tr <- if(html) {
              z <- cbind(Variable=paste0('$', varnames, '$'),
                         Transformation=paste0('$', TLi, '$'))
              as.character(htmlTable::htmlTable(z, caption='Pre-transformations',
                                                css.cell='min-width: 9em;',
                                                align='|l|l|',
                                                align.header='|c|c|',
                                                escape.html=FALSE))
            }
            else
              pretrans_tiny_table(varnames[i], TLi[i],
                                  if(prType() == 'typst') 'typst' else 'latex')
    }
    tex <- c(tex, tr)
  }
  if(file == '') return(tex)
  cat(tex, sep='\n', file=file, append=append)
}


## -----------------------------------------------------------------------
## NEW: latex_to_typst_math
## Targeted LaTeX-to-Typst converter for content ALREADY built (earlier
## in this function) using latexrms's own LaTeX-syntax from/to table --
## needed because that content is shared with the equation array (which
## goes through texmath and needs no translation) but is now ALSO fed to
## tinytable (which does NOT run cell content through texmath, confirmed
## directly by the "S_{0}(t)" braces-showing-literally bug on the
## survival table). Deliberately narrow: covers only the constructs
## latexrms's own from/to table is known to produce (\mathrm{}, \sqrt{},
## \log(, \min(, \max(, brace-based sub/superscripts), not a general
## LaTeX parser. \mathrm{} becomes upright(), not a quoted string, so
## nested math (e.g. a digit subscript already inside the \mathrm{}
## content) keeps working rather than being flattened to literal text --
## same reasoning as the mathvar design in markupSpecs$typst.
##
## This is the least-tested part of today's work -- regex-based reverse
## translation of syntax built elsewhere, not something compile-tested
## yet the way most of this effort has been. Worth a dedicated test
## before trusting it in place.
## -----------------------------------------------------------------------
latex_to_typst_math <- function(x) {
  x <- gsub('\\\\mathrm\\{([^}]*)\\}', 'upright(\\1)', x)
  x <- gsub('\\\\sqrt\\{([^}]*)\\}',   'sqrt(\\1)',    x)
  x <- gsub('\\\\log\\(', 'log(', x, fixed = TRUE)
  x <- gsub('\\\\min\\(', 'min(', x, fixed = TRUE)
  x <- gsub('\\\\max\\(', 'max(', x, fixed = TRUE)
  x <- gsub('_\\{([^}]*)\\}', '_(\\1)', x)
  x <- gsub('\\^\\{([^}]*)\\}', '^(\\1)', x)
  x
}


## -----------------------------------------------------------------------
## NEW: pretrans_tiny_table
## Pre-transformations table via tinytable, mirroring surv_tiny_table's
## structure. Both columns left-justified uniformly (no header-alignment
## override needed here, unlike the survival table). Title line uses the
## already-confirmed \textbf{} math-mode technique (same as the caption
## fix elsewhere) rather than depending on tinytable's own
## (unconfirmed) caption argument.
## -----------------------------------------------------------------------
pretrans_tiny_table <- function(vars, trans, output) {
  ## Local, non-exported helper -- see surv_tiny_table's identical
  ## definition for the full rationale.
  typstFence <- function(x) paste0("````{=typst}\n", x, "\n````\n")

  if(output == 'typst') {
    vars  <- latex_to_typst_math(vars)
    trans <- latex_to_typst_math(trans)
  }
  d <- data.frame(Variable       = paste0('$', vars,  '$'),
                  Transformation = paste0('$', trans, '$'))

  x <- tt(d, output = output)
  x <- style_tt(x, j = 1 : 2, align = 'l')

  if(output == 'typst')
    x <- theme_typst(x, figure = TRUE, align_figure = "c")
  else if(output == 'latex')
    ## Confirmed by direct compile test: theme_latex(placement=...)
    ## already inserts \centering automatically -- no separate
    ## centering mechanism needed for LaTeX.
    x <- theme_latex(x, placement = 'H')

  ## save_tt(x, output): documented public function, returns the
  ## rendered string directly. Confirmed reliable, unlike as.character()
  ## (real internal bug) and build_tt() (works but unexported).
  result <- save_tt(x, output)

  ## Same fence requirement as surv_tiny_table -- see its comment for
  ## the full explanation. typst needs the fence; latex stays bare.
  if(output == 'typst') result <- typstFence(result)

  title <- "$$\\textbf{Pre--Transformations}$$"
  c(title, result)
}
