latex.psm <-
  function(object,  title,
           file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
           append=FALSE, which=NULL, varnames, 
           columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='', ...)
{

  f <- object
  whichNot <- length(which)==0
  
  w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

  if(whichNot & !inline)
    {
      dist <- f$dist
      w <- c(w, paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",
                      survreg.auxinfo[[dist]]$latex(f$scale),
                      "{\\rm \\ \\ where} \\\\ \\]",sep=""))
    }
  atr <- f$Design

  if(whichNot) which <- 1:length(atr$name)
  if(missing(varnames)) varnames <- atr$name[atr$assume.code!=9]

  cat(w, sep=if(length(w)) "\n" else "", file=file, append=append)
  latexrms(f, file=file, append=TRUE, which=which,
           varnames=varnames, columns=columns, 
           before=before, after=after,
           prefix=if(whichNot)"X\\hat{\\beta}" else NULL, 
           inline=inline,pretrans=pretrans, digits=digits, size=size)
}


