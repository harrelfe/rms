cr.setup <- function(y)
{
  yname <- as.character(substitute(y))
  if(!is.factor(y)) y <- factor(y, exclude=NA)
  y       <- unclass(y)   # in case is.factor
  ylevels <- levels(y)
  kint    <- length(ylevels) - 1
  y       <- as.integer(y-1)
  
  reps   <- ifelse(is.na(y), 1, ifelse(y < kint-1, y+1, kint))
  subs   <- rep(1:length(y), reps)

  cuts   <- vector('list',kint+2)
  cuts[[1]] <- NA
  for(j in 0:kint) cuts[[j+2]] <- if(j < kint-1) 0:j else 0:(kint-1)

  cuts   <- unlist(cuts[ifelse(is.na(y),1,y+2)])
  y      <- rep(y, reps)
  Y      <- 1*(y==cuts)
  labels <- c('all', paste(yname,'>=',ylevels[2:kint],sep=''))
  cohort <- factor(cuts, levels=0:(kint-1), labels=labels)
  
  list(y=Y, cohort=cohort, subs=subs, reps=reps)
}

