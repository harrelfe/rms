latex.ols <-
  function(object, title,
           file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='', ...)
{
  f <- object
  
  w <- if(length(caption)) paste('\\begin{center} \\bf',
                                 caption,'\\end{center}')
  
  if(missing(which) & !inline)
    {
      Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
      
      w <- c(w, paste("\\[{\\rm E(",Y,
                      "}) = X\\beta, {\\rm \\ \\ where} \\\\ \\]", sep=""))
    }
  at <- f$Design
  
  if(missing(which)) which <- 1:length(at$name)
  
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  cat(w, file=file, sep=if(length(w)) "\n" else "", append=append)
  latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, after=after, prefix="X\\hat{\\beta}", inline=inline, 
           pretrans=pretrans, digits=digits, size=size)
}


