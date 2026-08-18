## -----------------------------------------------------------------------
## latex.lrm and latex.orm, with three changes in each, all marked below:
##
## CHANGE A (caption): \begin{center} \bf caption \end{center} -> a
## $$\textbf{caption}$$ math-mode construction. The old version sits
## outside any math delimiters, so it's plain text with embedded raw
## LaTeX commands -- Pandoc drops raw LaTeX silently for any output
## format other than Markdown/LaTeX/Org/ConTeXt (Typst isn't among
## them), meaning every caption would simply vanish under
## prType='typst'. Display math is center-aligned by default in both
## real LaTeX and Typst, so \textbf{} inside $$...$$ reproduces the same
## visual result while staying in the confirmed-working math-content
## pathway. \textbf{} in math mode confirmed by direct compile test:
## words render correctly and stay separated, though exact multi-space
## runs within the caption text are not preserved (harmless for this
## use -- captions don't depend on exact internal spacing).
##
## CHANGE B (array column spec): bare \begin{array} -> \begin{array}{lll}
## -- same bug already fixed in latexrms itself: texmath's array reader
## requires an explicit column-spec argument and errors without one.
## This array is 3 columns (alpha-subscript label, "=", coefficient),
## confirmed by the &=& structure in the row-building paste() calls
## below.
##
## CHANGE C (typo, unrelated to Typst): c(w,"\\end{array}",sep="") ->
## c(w,"\\end{array}") -- c() has no sep= parameter; the original was
## silently creating a spurious named element rather than doing
## anything useful. Noticed in passing, not a Typst-specific fix.
## -----------------------------------------------------------------------
latex.lrm <-
  function(object, title, 
           file='',
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "&", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           ...)
{
  f <- object
  md <- prType() %in% c('html', 'md', 'markdown')
  
  if(missing(which) & !inline)
    {
      Y <- paste("\\mathrm{", as.character(attr(f$terms,"formula")[2]),
                 "}", sep="")
      lev <- names(f$freq)
      nrp <- f$non.slopes
      
      w <- '$$'
      
      j <- if(lev[2]=="TRUE") "" else paste("=", lev[2], sep="")
      if(nrp==1) w <- paste(w, "P(", Y, j,
           ") = \\frac{1}{1+\\exp(-X\\beta)}", sep="")

      else
        w <- paste(w,"P(", Y, 
                   "\\geq j) = \\frac{1}{1+\\exp(-\\alpha_{j}-X\\beta)}",
                   sep="")

      w <- paste(w, ", \\mathrm{~~where}$$", sep="")

      if(length(caption)) {
        if(md) w <- c(paste('<div align=center><strong>', caption,
                             '</strong></div>'), w)
        else
          w <- c(paste0("$$\\textbf{", caption, "}$$"), w)   ## CHANGE A
        }
      
      if(nrp > 1) {
        w <- c(w,"$$\\begin{array}{lll}")   ## CHANGE B, now $$-wrapped
        cof <- format(f$coef[1:nrp])
        for(i in 1:nrp)
          w <- c(w, paste("\\hat{\\alpha}_{\\rm ",
                          lev[i+1],"} &=&",cof[i],"\\\\",sep=""))
        w <- c(w,"\\end{array}$$")   ## CHANGE C, now $$-wrapped
      }
    }
  else w <- NULL
  
  if(missing(which) | missing(varnames)) at <- f$Design

  if(missing(which))    which    <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]

  z <- latexrms(f, file='', which=which, varnames=varnames, 
           columns=columns, 
           before=before, after=after, prefix="X\\hat{\\beta}",
           inline=inline, pretrans=pretrans, digits=digits,
           size=size)
  if(inline) return(z)
  w <- c(w, z)
  if(file != '' || prType() == 'plain') {
    cat(w, file=file, append=append, sep='\n')
    return(invisible())
    }
  rendHTML(w, html=FALSE)
}


latex.orm <-
  function(object, title, 
           file='',
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "&", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           intercepts=nrp < 10, ...)
{
  f <- object

  md <- prType() %in% c('html', 'md', 'markdown')

  
  if(missing(which) & !inline)
    {
      Y <- paste0("\\mathrm{", f$yname, "}")
      lev <- f$yunique
      nrp <- f$non.slopes

      z <- '\\alpha_{y} + X\\beta'
      zm <- '- \\alpha_{y} - X\\beta'
      dist <-
        switch(f$family,
               logistic = paste('\\frac{1}{1+\\exp(', zm, ')}', sep=''),
               probit   = paste('\\Phi(', z, ')', sep=''),
               cauchit  = paste('\\frac{1}{\\pi}\\tan^{-1}(', z,
                 ') + \\frac{1}{2}', sep=''),
               loglog   = paste('\\exp(-\\exp(', zm, '))', sep=''),
               cloglog  = paste('1 - \\exp(-\\exp(', z, ')', sep=''))
                     
      w <- '$$'
      
      w <- paste(w, "P(", Y, 
                   "\\geq y | X) = ", dist, sep='')

      w <- paste(w, "\\mathrm{~~where}$$", sep="")

      if(length(caption)) {
        if(md) w <- c(paste('<div align=center><strong>', caption,
                             '</strong></div>'), w)
        else
          w <- c(paste0("$$\\textbf{", caption, "}$$"), w)   ## CHANGE A
        }
      
      if(intercepts) {
        nl <- as.numeric(lev)
        if(!any(is.na(nl))) lev <- format(nl, digits=digits)
          w <- c(w,"$$\\begin{array}{lll}")   ## CHANGE B, now $$-wrapped
          cof <- format(f$coef[1:nrp], digits=digits)
          for(i in 1:nrp)
            w <- c(w, paste("\\hat{\\alpha}_{\\mathrm{",
                            lev[i+1], "}} &=&", cof[i], "\\\\", sep=""))
          w <- c(w, "\\end{array}$$")   ## CHANGE C, now $$-wrapped
      }
    }
  else w <- NULL
  if(missing(which) | missing(varnames)) at <- f$Design

  if(missing(which)) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  z <- 
  latexrms(f, file='', append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, after=after, prefix="X\\hat{\\beta}",
           inline=inline, pretrans=pretrans, digits=digits,
           size=size)
  if(inline) return(z)
  w <- c(w, z)
  if(file == '' && prType() != 'plain') return(rendHTML(w, html=FALSE))
  cat(w, file=file, append=append, sep='\n')
  invisible()
}
