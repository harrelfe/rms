ie.setup <- function(failure.time, event, ie.time, break.ties=FALSE)
  {
    s <- !is.na(ie.time)
    if(all(s)) warning('every subject had an intervening event')
    if(!any(s)) stop('no intervening events')
    if(any(ie.time[s] > failure.time[s])) 
      stop('an ie.time was after a failure.time')
    if(break.ties)
      {
        mindif <-
          min(diff(sort(unique(failure.time[!is.na(failure.time)]))))
        ## 8Nov01 Thanks: Josh Betcher
        k <- s & (ie.time==failure.time)
        if(sum(k)==0) warning('break.times=T but no ties found')
        ie.time[k] <- ie.time[k] - runif(sum(k),0,mindif)
      }
 
    if(any(ie.time[s]==failure.time[s])) 
      stop('ie.times not allowed to equal failure.times')

    n <- length(failure.time)
    reps <- ifelse(is.na(ie.time), 1, 2)
    subs <- rep(1:n, reps)

    start <- end <- ev <- ie.status <- vector('list', n)
    start[]     <- 0
    end[]       <- failure.time
    ev[]        <- event
    ie.status[] <- 0
    for(i in seq(along=s)[s])
      {
        start[[i]]  <- c(0, ie.time[i])
        end[[i]]    <- c(ie.time[i], failure.time[i])
        ev[[i]]     <- c(0, event[i])
        ie.status[[i]] <- c(0, 1)
      }

    start     <- unlist(start)
    end       <- unlist(end)
    ev        <- unlist(ev)
    ie.status <- unlist(ie.status)
    
    u <- units(failure.time)
    units(end) <- if(u=='')'Day' else u
    
    s <- !is.na(start+end) & (end <= start)
    if(any(s))
      {
        cat('stop time <= start time:\n')
        print(cbind(start=start[s], end=end[s]))
        stop()
      }

    S <- Surv(start, end, ev)

    list(S=S, ie.status=ie.status, subs=subs, reps=reps)
  }
