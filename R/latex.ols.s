latex.ols <-
  function(object, title,
           file='',
           append=FALSE, which, varnames, columns=65, inline=FALSE, 
           before=if(inline)"" else "& &", after="",
           pretrans=TRUE, caption=NULL, digits=.Options$digits, size='',
           ...)
{
  f <- object

  md <- prType() %in% c('html', 'md', 'markdown')
  
  w <- if(length(caption)) {
         if(md) paste('<div align=center><strong>', caption,
                      '</strong></div>', sep='')
         else
           paste('\\begin{center} \\bf',
                 caption,'\\end{center}')
         }
  
  if(missing(which) & !inline)
    {
      Y <- paste("{\\rm ",as.character(attr(f$terms,"formula")[2]),"}",sep="")
      
      w <- c(w, paste("\\[{\\rm E(",Y,
                      "}) = X\\beta, {\\rm \\ \\ where} \\\\ \\]", sep=""))
    }
  at <- f$Design
  
  if(missing(which)) which <- 1:length(at$name)
  
  if(missing(varnames)) varnames <- at$name[at$assume.code!=9]
  if(! md)
    cat(w, file=file, sep=if(length(w)) "\n" else "", append=append)
  z <- latexrms(f, file=file, append=TRUE, which=which, varnames=varnames, 
                columns=columns, 
                before=before, after=after, prefix="X\\hat{\\beta}",
                inline=inline, 
                pretrans=pretrans, digits=digits, size=size)
  if(md) htmltools::HTML(c(paste0(w, '\n'), as.character(z)))
  else z
}


