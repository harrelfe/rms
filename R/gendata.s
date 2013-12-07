gendata <- function(fit, ..., nobs, viewvals=FALSE, expand=TRUE, factors)
  {
    at <- fit$Design
    nam <- at$name[at$assume!="interaction"]

    if(!missing(nobs) && !is.logical(nobs))
      {
        df <- predictrms(fit, type="adjto.data.frame")
        df[1:nobs,] <- df
        cat("Edit the list of variables you would like to vary.\nBlank out variables to set to reference values.\n")
        nam.sub <- de(nam)[[1]]
        nam.sub <- nam.sub[!is.na(nam.sub)]
        if(!all(nam.sub %in% nam)) stop("misspelled a variable name")
        df.sub <- as.data.frame(df[,nam.sub])
        cat("Edit the predictor settings to use.\n")
        if(viewvals && 
           length(val <- Getlim(at, allow.null=TRUE,
                                need.all=FALSE)$values[nam.sub]))
          {
            cat("A window is being opened to list the valid values of discrete variables.\n")
            sink(tf <- tempfile())
            print.datadist(list(values=val))
            sink()
            file.show(tf)
          }
        for(i in 1:length(df.sub))
          if(is.factor(df.sub[[i]])) df.sub[[i]] <- as.character(df.sub[[i]])
        df.sub <- as.data.frame(de(df.sub))
        df[nam.sub] <- df.sub
        return(structure(df, names.subset=nam.sub))
      }

    factors <- if(missing(factors)) rmsArgs(substitute(list(...))) else factors
    fnam <- names(factors)
    nf <- length(factors)
    if(nf==0) return(predictrms(fit, type="adjto.data.frame"))
    which <- charmatch(fnam, nam, 0)
    if(any(which==0)) stop(paste("factor(s) not in design:",
             paste(names(factors)[which==0],collapse=" ")))
    settings <- if(nf<length(nam)) predictrms(fit, type="adjto.data.frame")
     else list()
    settings <- unclass(settings)
    if(nf>0) for(i in 1:nf) settings[[fnam[i]]] <- factors[[i]]
    if(nf==0) return(as.data.frame(settings))
    if(expand) expand.grid(settings) else as.data.frame(settings)
  }
