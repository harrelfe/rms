latex.lrm <- function(object, title, 
   file=paste(first.word(deparse(substitute(object))),".tex",sep=""),
   append=FALSE, which, varnames, columns=65, inline=FALSE, 
   before=if(inline)"" else "& &",pretrans=TRUE,caption=NULL,...)
{
  f <- object
  
  if(missing(which) & !inline)
    {
      Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
      lev <- names(f$freq)
      nrp <- f$non.slopes
      
      w <- '\\['
      
      j <- if(lev[2]=="TRUE") "" else paste("=",lev[2],sep="")
      if(nrp==1) w <- paste(w,"{\\rm Prob}\\{",Y, j,
           "\\} = \\frac{1}{1+\\exp(-X\\beta)}", sep="")

      else
        w <- paste(w,"{\\rm Prob}\\{", Y, 
                   "\\geq j\\} = \\frac{1}{1+\\exp(-\\alpha_{j}-X\\beta)}",
                   sep="")

      w <- paste(w, ", {\\rm \\ \\ where} \\\\ \\]", sep="")

      if(length(caption)) w <- c(paste('\\begin{center} \\bf',caption,
                                       '\\end{center}'), w)
      
      if(nrp>1)
        {
          w <- c(w,"\\begin{eqnarray*}")
          cof <- format(f$coef[1:nrp])
          for(i in 1:nrp)
            w <- c(w, paste("\\hat{\\alpha}_{\\rm ",
                            lev[i+1],"} &=&",cof[i],"\\\\",sep=""))
          w <- c(w,"\\end{eqnarray*}",sep="")
        }
							}
  else w <- NULL
  if(missing(which) | missing(varnames)) at <- f$Design

  if(missing(which)) which <- 1:length(at$name)
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  cat(w, file=file, append=append, sep=if(length(w))"\n" else "")
  latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
           columns=columns, 
           before=before, prefix="X\\hat{\\beta}", inline=inline,
           pretrans=pretrans)
}


