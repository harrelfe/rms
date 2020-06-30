which.influence <- function(fit, cutoff=.2)
{
  cox <- inherits(fit,"cph")
  
  stats <- resid(fit, "dfbetas")
  rnam  <- which(! is.na(stats[,1]))
  stats <- stats[rnam,, drop=FALSE]   ##delete rows added back due to NAs
  d <- dimnames(stats)[[1]]
  if(length(d)) rnam <- d
  
  at <- fit$Design
  
  w      <- list()
  namw   <- NULL
  k      <- 0
  oldopt <- options('warn')
  options(warn=-1)
  on.exit(options(oldopt))
  
  if(! cox)
    {
      ww <- rnam[abs(stats[, 1]) >= cutoff]
      if(length(ww))
        {
          k      <- k + 1
          w[[k]] <- ww
          namw   <- "Intercept"
        }
    }
  
  Assign <- fit$assign
  nrp    <- num.intercepts(fit)
  assadj <- if(nrp > 1) nrp - 1 else 0
  nm     <- names(Assign)[1]
  if(nm=="Intercept" | nm=="(Intercept)") Assign[[1]] <- NULL
  ##remove and re-number

  j <- 0
  for(i in (1 : length(at$name))[at$assume.code != 8])
    {
      j <- j + 1
      as <- Assign[[j]] - assadj
      if(length(as) == 1) ww <- rnam[abs(stats[, as]) >= cutoff]
	  else
        {
          z <- rep(FALSE, length(rnam))
          for(r in as)
            z <- z | abs(stats[, r]) >= cutoff
            ww <- rnam[z]
        }
      if(length(ww))
        {
          k      <- k+1
          w[[k]] <- ww
          namw   <- c(namw, at$name[i])
        }
     }
  if(length(w)) names(w) <- namw
  
  w
}


##show.influence was written by:
##Jens Oehlschlaegel-Akiyoshi
##oehl@psyres-stuttgart.de
##Center for Psychotherapy Research
##Christian-Belser-Strasse 79a
##D-70597 Stuttgart Germany

show.influence <- function(object, dframe, report=NULL, sig=NULL, id=NULL)
{
  who <- unlist(object)
  nam <- names(object)
  ## In future parse out interaction components in case main effects
  ## not already selected
  ia <- grep('\\*', nam)              # remove interactions
  if(length(ia)) nam <- nam[-ia]
  nam  <- nam[nam %nin% 'Intercept']  # remove Intercept
  rnam <- dimnames(dframe)[[1]]
  if(! length(rnam)) rnam <- 1:nrow(dframe)
  if (length(report)) col <- c(nam,
                               dimnames(dframe[,report,drop=FALSE])[[2]] )
  else col <- nam
  row <- rnam %in% who
  if(any(col %nin% names(dframe)))
    stop(paste('needed variables not in dframe:',
               paste(col[col %nin% names(dframe)],collapse=' ')))
  dframe <- dframe[row,col,drop=FALSE]
  rnam   <- rnam[row]
  Count  <- table(who)
  Count  <- as.vector(Count[match(rnam,names(Count))])
  for (i in 1 : length(nam))
    {
      ni  <- nam[i]
      val <- dframe[,ni]
      if (length(sig) && is.numeric(val)) val <- signif(val, sig) else
      val <- format(val)
      dframe[,ni] <- paste(ifelse(rnam %in% object[[ni]],"*",""), val, sep  = "")
      ## In future change i to also find any object containing the
      ## variable (e.g., interaction)   was object[[i]] dframe[,i] 24Nov00
    }
  if (length(sig) && length(report))
	for (i in (length(nam) + 1) : dim(dframe)[2])
      if(is.numeric(dframe[, i]))
        dframe[,i] <- signif(dframe[, i], sig)
  dframe <- data.frame(Count,dframe)
  if(length(id)) row.names(dframe) <- id[as.numeric(row.names(dframe))]
  print(dframe, quote=FALSE)
  invisible(dframe)
}

