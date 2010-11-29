latex.pphsm <- function(object, title,
    file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
    append=FALSE, which=NULL, varnames, 
    columns=65, inline=FALSE, 
    before=if(inline)"" else "& &",pretrans=TRUE, caption=NULL, ...)
{
  whichThere <- length(which)
  w <- if(length(caption)) paste('\\begin{center} \\bf',caption,'\\end{center}')

  sc <- exp(object$parms)
  at <- object$Design

  if(!whichThere & !inline)
    {
      dist <- paste("\\exp\\{-t^{",format(1/sc),"} \\exp(X\\hat{\\beta})\\}")
      w <- c(w,paste("\\[{\\rm Prob}\\{T\\geq t\\} = ",dist,
                     "{\\rm \\ \\ where} \\\\ \\]",sep=""))
    }				
  if(!whichThere) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  cat(w, file=file, sep=if(length(w))"\n" else "", append=append)
  latexrms(object, file=file, append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, prefix=if(!whichThere)"X\\hat{\\beta}" else NULL, 
           inline=inline,pretrans=pretrans)
}


