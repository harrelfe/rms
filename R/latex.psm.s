latex.psm <-
  function(object,  title,
           file='',
           append=FALSE, which=NULL, varnames, 
           columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           ...)
{
  md <- prType() %in% c('html', 'md', 'markdown')
  
  f <- object
  whichNot <- length(which)==0
  
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>', sep='')
         else
           paste('\\begin{center} \\bf',caption,'\\end{center}')
         }

  if(whichNot & !inline)
    {
      dist <- f$dist
      w <- c(w, paste("$$\\Pr(T\\geq t) = ",
                      survreg.auxinfo[[dist]]$latex(f$scale),
                      "~\\mathrm{where}$$",sep=""))
    }
  atr <- f$Design

  if(whichNot) which <- 1:length(atr$name)
  if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]

  if(file != '') cat(w, sep=if(length(w)) "\n" else "",
                     file=file, append=append)
  ltx <- latexrms(f, append=TRUE, which=which,
                  varnames=varnames, columns=columns, 
                  before=before, after=after,
                  prefix=if(whichNot)"X\\hat{\\beta}" else NULL, 
                  inline=inline,pretrans=pretrans, digits=digits,
                  size=size)
  if(inline) return(ltx)
  z <- c(w, ltx)
  if(file == '' && prType() != 'plain') return(rendHTML(z, html=FALSE))
  cat(z, file=file, append=append, sep='\n')
  invisible()
}
